#!/usr/bin/env python3

"""Filter per-sample promoter DM calls; annotate POG calls with Ploidetect.

Requires the public promoter_methylation_nanopore directory layout.
Original thresholds and order are retained. Does not use Andrew's analysis.
"""
import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

LOG = logging.getLogger(__name__)
COMPARISON_SUFFIX = '_dmr_comparison_to_cohort.tsv.gz'
OUT_SUFFIX = '_variants_filtered.tsv.gz'
CN_COLS = ('zygosity', 'LOH', 'A', 'B', 'segment_depth')


def read_tsv(path, **kwargs):
    if not path.is_file():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep='\t', compression='infer', **kwargs)


def methylation_positions(path):
    """Keep every CpG row, including duplicate positions (historical behavior)."""
    if not path.is_file():
        raise FileNotFoundError(path)
    # Methylation frequency files are BED-like, without a header.
    data = read_tsv(path, header=None, usecols=[0, 1], names=['chrom', 'start'])
    data['start'] = pd.to_numeric(data['start'], errors='raise').astype('int64')
    return {chrom: np.sort(group['start'].to_numpy())
            for chrom, group in data.groupby('chrom', sort=False)}


def count_positions(index, chrom, start, end):
    positions = index.get(chrom)
    if positions is None:
        raise ValueError(f'Chromosome {chrom} missing from CpG frequency file')
    # Inclusive start/end as in original script.
    return int(np.searchsorted(positions, end, side='right') -
               np.searchsorted(positions, start, side='left'))


def frequency_path(root, sample_type, haplotype, sample):
    folder = root / 'methylation_frequency' / sample_type / haplotype
    matches = sorted(folder.glob(f'{sample}_{haplotype}_methyl_freq_intersect_w_filtered_dss_regions.tsv*'))
    if len(matches) != 1:
        raise FileNotFoundError(f'Expected exactly one {haplotype} frequency file for {sample} in {folder}; found {len(matches)}')
    return matches[0]


def annotate_cn(frame, path):
    cn = read_tsv(path, low_memory=False)
    gene_col = 'gene' if 'gene' in cn else 'gene_name' if 'gene_name' in cn else None
    if gene_col is None:
        raise ValueError(f'Ploidetect file missing gene/gene_name: {path}')
    cn = cn.drop_duplicates(subset=[gene_col], keep='first').set_index(gene_col)
    for col in CN_COLS:
        frame[col] = frame['dmr_hugo'].map(cn[col]) if col in cn else np.nan
    return frame


def load_comparison(path):
    df = read_tsv(path)
    if 'dmr_region' not in df:
        if 'Unnamed: 0' in df:
            df = df.rename(columns={'Unnamed: 0': 'dmr_region'})
        elif df.columns[0].startswith('chr') and df.columns[0] not in ('chr', 'chrom'):
            df = df.rename(columns={df.columns[0]: 'dmr_region'})
        else:
            raise ValueError(f'Cannot identify DMR region identifier column: {path}')
    required = {'dmr_hugo', 'nCpG_dmr_file', 'nCpG_methyl_freq_file', 'perc_CpG_outlier'}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f'{path} missing columns {sorted(missing)}')
    return df


def process_sample(path, sample, sample_type, root, outdir, reftss, imprinted):
    df = load_comparison(path)
    initial = len(df)
    # Identifier convention: chr_start_end_gene (gene may contain underscores).
    parts = df['dmr_region'].astype(str).str.split('_', n=3, expand=True)
    if parts.shape[1] < 4:
        raise ValueError(f'Invalid DMR identifiers in {path}')
    df['region_start'] = pd.to_numeric(parts[1], errors='raise').astype('int64')
    df['region_end'] = pd.to_numeric(parts[2], errors='raise').astype('int64')
    df['region_length'] = (df['region_end'] - df['region_start']).abs()
    df['chrX_chrY'] = np.where(parts[0].isin(['chrX', 'chrY']), 'yes', 'no')
    df['CpG_outlier'] = np.where(pd.to_numeric(df['perc_CpG_outlier']) >= .75, 'yes', 'no')
    denom = pd.to_numeric(df['nCpG_dmr_file'])
    numerator = pd.to_numeric(df['nCpG_methyl_freq_file'])
    df['perc_CpG_found'] = np.where(denom > 0, numerator * 100 / denom, np.nan)
    df = df.loc[df['chrX_chrY'].eq('no')].copy()
    n_autosomal = len(df)
    df = df.loc[df['CpG_outlier'].eq('yes')].copy()
    n_outlier = len(df)
    df['found_more_than_75perc'] = np.where(df['perc_CpG_found'] >= 75, 'yes', 'no')
    df = df.loc[df['found_more_than_75perc'].eq('yes')].copy()
    n_found = len(df)

    hp1 = methylation_positions(frequency_path(root, sample_type, 'HP1', sample)) if len(df) else {}
    hp2 = methylation_positions(frequency_path(root, sample_type, 'HP2', sample)) if len(df) else {}
    promoter_counts = []
    for gene in df['dmr_hugo']:
        if gene not in reftss.index:
            promoter_counts.append(np.nan)
            continue
        tss = reftss.loc[gene]
        if isinstance(tss, pd.DataFrame):
            raise ValueError(f'Duplicate refTSS symbol {gene}; original expected unique symbols')
        tss_start = int(tss['start'])
        start, end = ((tss_start - 1500, tss_start + 1000)
                      if tss['strand'] == '+' else
                      (tss_start - 1000, tss_start + 1500))
        promoter_counts.append(max(count_positions(hp1, tss['chrom'], start, end),
                                   count_positions(hp2, tss['chrom'], start, end)))
    df['#ofCpG_in_promoter'] = promoter_counts
    df['30perc_of_promoter_CpG'] = df['#ofCpG_in_promoter'] * .30
    df['cover_more_than_30perc_promoter_CpG'] = np.where(
        df['#ofCpG_in_promoter'].notna() &
        (pd.to_numeric(df['nCpG_methyl_freq_file']) >= df['30perc_of_promoter_CpG']),
        'yes', 'no')
    df = df.loc[df['cover_more_than_30perc_promoter_CpG'].eq('yes')].copy()
    # Historical hypermethylation count is HERE, before imprinted-gene exclusion.
    pre_imprinted_count = len(df)
    df['imprinted'] = np.where(df['dmr_hugo'].isin(imprinted), 'yes', 'no')
    if sample_type == 'POG':
        cn_path = root / 'copy_number' / 'POG' / f'{sample}_cna_genes.bed.gz'
        df = annotate_cn(df, cn_path)
    df = df.loc[df['imprinted'].eq('no')].copy()
    outdir.mkdir(parents=True, exist_ok=True)
    df.to_csv(outdir / f'{sample}{OUT_SUFFIX}', sep='\t', index=False, compression='gzip')
    return dict(tumour_original_source=sample, input_rows=initial,
                autosomal_rows=n_autosomal, outlier_rows=n_outlier,
                cpg_representation_rows=n_found, total_count=pre_imprinted_count,
                final_filtered_count=len(df))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample-type', choices=['POG', 'BLOOD'], required=True)
    p.add_argument('--data-dir', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--sample', help='Optional single sample; full-cohort run needed for hypermethylation threshold')
    args = p.parse_args()
    root = args.data_dir
    reftss = read_tsv(root / 'annotations/refTSS/refTSS_hg38_single_best_transcript.bed.gz').set_index('symbol')
    imprinted = set(read_tsv(root / 'annotations/imprinted_genes/hg38_mapped_imprinted_genes_geneimprint_db.tsv.gz')['hugo'].dropna())
    folder = root / 'dmr_comparison_to_cohort' / args.sample_type
    paths = sorted(folder.glob(f'*{COMPARISON_SUFFIX}'))
    if args.sample:
        paths = [x for x in paths if x.name == f'{args.sample}{COMPARISON_SUFFIX}']
    if not paths:
        p.error(f'No matching comparison files in {folder}')
    records = []
    for path in paths:
        sample = path.name[:-len(COMPARISON_SUFFIX)]
        LOG.info('Filtering %s', sample)
        records.append(process_sample(path, sample, args.sample_type, root,
                                      args.output_dir / 'per_sample', reftss, imprinted))
    counts = pd.DataFrame(records).sort_values('tumour_original_source')
    counts_dir = args.output_dir / 'counts'
    counts_dir.mkdir(parents=True, exist_ok=True)
    counts.to_csv(counts_dir / 'summary_counts_per_sample.tsv', sep='\t', index=False)
    if args.sample_type == 'POG':
        counts[['tumour_original_source', 'total_count']].to_csv(
            counts_dir / 'summary_counts_per_pog.tsv', sep='\t', index=False)
        if not args.sample:
            mean = counts['total_count'].mean()
            sd = counts['total_count'].std(ddof=0)
            threshold = mean + 2 * sd
            excluded = counts.loc[counts['total_count'] >= threshold].copy()
            excluded.to_csv(counts_dir / 'hypermethylated_samples.tsv', sep='\t', index=False)
            pd.DataFrame([{'n_samples': len(counts), 'mean': mean, 'std_ddof0': sd,
                           'threshold_mean_plus_2std': threshold}]).to_csv(
                counts_dir / 'hypermethylation_threshold.tsv', sep='\t', index=False)
            LOG.info('POG hypermethylation threshold %.4f; excluded %d samples', threshold, len(excluded))
    LOG.info('Processed %d samples', len(counts))


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    main()
