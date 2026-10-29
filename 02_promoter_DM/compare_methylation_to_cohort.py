#!/usr/bin/env python3

"""Compare per-sample haplotype CpG methylation differences with cohort references."""
import argparse
from pathlib import Path
import pandas as pd
import numpy as np

SUFFIX = '_methyl_freq_intersect_w_filtered_dss_regions.tsv.gz'
COLUMNS = ['chr', 'start', 'end', 'all', 'methylated', 'freq']
OUT_COLUMNS = ['dmr_hugo', 'meanMethy1_dmr_file', 'meanMethy2_dmr_file', 'nCpG_dmr_file', 'nCpG_methyl_freq_file', 'nCpG_outlier', 'perc_CpG_outlier']

def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample-type', choices=['POG', 'BLOOD'], required=True)
    p.add_argument('--data-dir', type=Path, required=True, help='Parent directory containing dss/ and methylation_frequency/')
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--overwrite', action='store_true')
    return p.parse_args()

def sample_files(directory, haplotype):
    suffix = f'_{haplotype}{SUFFIX}'
    return {p.name[:-len(suffix)]: p for p in (directory / haplotype).glob(f'*{suffix}')}

def read_cpg(path):
    df = pd.read_csv(path, sep='\t', header=None, names=COLUMNS, compression='infer')
    df['start'] = pd.to_numeric(df['start'], errors='raise')
    df['end'] = pd.to_numeric(df['end'], errors='raise')
    df['freq'] = pd.to_numeric(df['freq'], errors='raise')
    if df.duplicated(['chr', 'start', 'end']).any():
        raise ValueError(f'Duplicate CpG coordinates in {path}; resolve before comparison')
    return df[['chr', 'start', 'end', 'freq']]

def main():
    args = parse_args()
    base = args.data_dir / "methylation_frequency" / args.sample_type
    dmr_file = args.data_dir / "dss" / args.sample_type / "filtered_dss_refTSS_hg38_single_best_transcript_intersect_50percent_overlap.tsv.gz"
    hp1_files = sample_files(base, 'HP1')
    hp2_files = sample_files(base, 'HP2')
    samples = set(hp1_files) & set(hp2_files)
    dmr = pd.read_csv(dmr_file, sep='\t')
    needed = {'POGID', 'chr', 'start', 'end', 'nCG', 'meanMethy1', 'meanMethy2', 'symbol'}
    if not needed.issubset(dmr.columns):
        raise ValueError(f'DMR table missing columns: {sorted(needed - set(dmr.columns))}')
    dmr['start'] = pd.to_numeric(dmr['start'])
    dmr['end'] = pd.to_numeric(dmr['end'])
    args.output_dir.mkdir(parents=True, exist_ok=True)
    reference_cache = {}
    print(f'{args.sample_type}: {len(samples)} paired samples; {dmr.POGID.nunique()} samples in DMR table')
    for sample, group in dmr.groupby('POGID', sort=True):
        if sample not in samples:
            print(f'WARNING: skipping {sample}: HP1 or HP2 file missing')
            continue
        output = args.output_dir / f'{sample}_dmr_comparison_to_cohort.tsv'
        if output.exists() and not args.overwrite:
            print(f'Skipping existing: {output}')
            continue
        hp1 = read_cpg(hp1_files[sample]).rename(columns={'freq': 'HP1_freq'})
        hp2 = read_cpg(hp2_files[sample]).rename(columns={'freq': 'HP2_freq'})
        paired = hp1.merge(hp2, on=['chr', 'start', 'end'], how='inner', validate='one_to_one')
        paired['diff'] = (paired.HP1_freq - paired.HP2_freq).abs()
        results = {}
        for chrom, chrom_group in group.groupby('chr', sort=False):
            if chrom not in reference_cache:
                ref_path = base / 'haplotype_difference_statistics' / f'{chrom}_mean_std_methyl_freq_diff_cpgs.tsv.gz'
                if not ref_path.exists():
                    raise FileNotFoundError(ref_path)
                ref = pd.read_csv(ref_path, sep='\t', compression='infer')
                ref = ref.rename(columns={ref.columns[0]: 'CpG'})
                required = {'CpG', 'mean_methyl_freq', 'std_methyl_freq'}
                if not required.issubset(ref.columns):
                    raise ValueError(f'Missing columns in {ref_path}')
                if ref.CpG.duplicated().any():
                    raise ValueError(f'Duplicate reference CpGs in {ref_path}')
                reference_cache[chrom] = ref
            chrom_cpg = paired.loc[paired.chr == chrom].copy()
            chrom_cpg['CpG'] = (chrom_cpg.chr.astype(str) + '_' + chrom_cpg.start.astype(str) + '_' + chrom_cpg.end.astype(str))
            chrom_cpg = chrom_cpg.merge(reference_cache[chrom], on='CpG', how='inner', validate='one_to_one')
            chrom_cpg['outlier'] = chrom_cpg['diff'] > 2 * chrom_cpg['std_methyl_freq']
            for row in chrom_group.itertuples(index=False):
                selected = chrom_cpg.loc[chrom_cpg.start.between(row.start, row.end, inclusive='both')]
                count = len(selected)
                n_outlier = int(selected.outlier.sum())
                key = f'{row.chr}_{row.start}_{row.end}_{row.symbol}'
                if key not in results:
                    results[key] = [row.symbol, row.meanMethy1, row.meanMethy2, row.nCG,
                                    count, n_outlier, n_outlier / count if count else 0]
        out = pd.DataFrame.from_dict(results, orient='index', columns=OUT_COLUMNS)
        out.to_csv(output, sep='\t', index=True)
        print(f'{sample}: {len(out)} DMR-promoter records -> {output}')

if __name__ == '__main__':
    main()
