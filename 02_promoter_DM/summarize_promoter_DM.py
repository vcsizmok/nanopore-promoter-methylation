#!/usr/bin/env python3

"""Summarize per-sample promoter DM and apply normal, lymphoid and HGNC filters."""
import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd

SUFFIX = '_variants_filtered.tsv.gz'

# Original study exclusion set; fixed for exact cohort reproduction.
HYPERMETHYLATED_POG_SAMPLES = {
    'POG1297-OCT-1', 'POG569-FF-1', 'POG1780-OCT-1',
    'POG279-OCT-1', 'POG1771-OCT-1', 'POG707-OCT-1',
    'POG1772-OCT-1', 'POG785-OCT-1', 'POG196-OCT-1',
    'POG688-OCT-1', 'POG1756-OCT-1', 'POG1022-OCT-1',
    'POG1779-OCT-1', 'POG1863-OCT-1', 'POG1098-OCT-1',
}

# Original post-imprinted filtering counts, for an optional audit.
ORIGINAL_HYPERMETHYLATED_COUNTS = {
    'POG1297-OCT-1': 677, 'POG569-FF-1': 548,
    'POG1780-OCT-1': 513, 'POG279-OCT-1': 441,
    'POG1771-OCT-1': 420, 'POG707-OCT-1': 385,
    'POG1772-OCT-1': 354, 'POG785-OCT-1': 354,
    'POG196-OCT-1': 328, 'POG688-OCT-1': 312,
    'POG1756-OCT-1': 306, 'POG1022-OCT-1': 295,
    'POG1779-OCT-1': 292, 'POG1863-OCT-1': 288,
    "POG1098-OCT-1": 291,
}

BASE_MAP = {
    'C11orf80': 'TOP6BL', 'C12orf40': 'REDIC1', 'C2orf16': 'SPATA31H1',
    'C7orf31': 'CYCS', 'FAM166A': 'CIMIP2A', 'FAM166C': 'CIMIP2C',
    'FAM183BP': 'CFAP144P1', 'KHDRBS2-OT1': 'KHDRBS2',
    'LINC00461': 'MIR9-2HG', 'METTL7A': 'TMT1A', 'ODF3': 'CIMAP1A',
    'PIFO': 'CIMAP3', 'RBM8B': 'RMB8AP1', 'RNF165': 'ARK2C',
    'SOGA3': 'MTCL3', 'TEX37': 'SPMIP9', 'TEX45': 'SAXO5',
    'TEX49': 'SPMIP11', 'THEGL': 'SPMAP2L', 'TTC26': 'IFT56',
    'C4orf48': 'NICOL1', 'ODF3B': 'CIMAP1B',
}


def read(path, **kwargs):
    if not path.is_file():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep='\t', compression='infer', **kwargs)


def save(df, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep='\t', index=False)
    print(f'{path}: {len(df)} genes')


def summarize(paths, sample_type, excluded):
    """Count only HET/HOM POG observations; audit unclassified observations."""
    unclassified = []
    if sample_type == 'POG':
        events = defaultdict(lambda: defaultdict(list))
    else:
        events = defaultdict(list)
    for path in paths:
        sample = path.name[:-len(SUFFIX)]
        if sample in excluded:
            continue
        df = read(path)
        if 'dmr_hugo' not in df:
            raise ValueError(f'{path} missing dmr_hugo')
        if sample_type == 'POG' and 'zygosity' not in df:
            print(f'WARNING: {path.name} has no zygosity column; auditing all observations')
            df['zygosity'] = pd.NA
        for _, row in df.iterrows():
            gene = row['dmr_hugo']
            if pd.isna(gene) or any(x in str(gene) for x in ('-AS1', '-AS2', '-AS3')):
                continue
            if sample_type == 'POG':
                z = row['zygosity']
                if pd.isna(z) or str(z).strip().lower() in ('', 'nan', 'none', 'na'):
                    unclassified.append({'sample': sample, 'gene': gene,
                                         'zygosity': '', 'reason': 'missing_zygosity'})
                    continue
                z = str(z).strip()
                if z not in ('HET', 'HOM'):
                    unclassified.append({'sample': sample, 'gene': gene,
                                         'zygosity': z, 'reason': 'unexpected_zygosity'})
                    continue
                events[gene][z].append(sample)
            else:
                events[gene].append(sample)
    rows = []
    for gene, value in events.items():
        if sample_type == 'POG':
            het, hom = value.get('HET', []), value.get('HOM', [])
            rows.append({'gene': gene, 'total_count': len(het) + len(hom),
                         'pogs_w_DMs': ', '.join(het + hom),
                         'aDM_count': len(het), 'pogs_w_aDMs': ', '.join(het),
                         'hAM_count': len(hom), 'pogs_w_hAMs': ', '.join(hom)})
        else:
            rows.append({'gene': gene, 'total_count': len(value),
                         'pogs_w_DMs': ', '.join(value)})
    columns = (['gene', 'total_count', 'pogs_w_DMs', 'aDM_count', 'pogs_w_aDMs',
                'hAM_count', 'pogs_w_hAMs'] if sample_type == 'POG'
               else ['gene', 'total_count', 'pogs_w_DMs'])
    return pd.DataFrame(rows, columns=columns), pd.DataFrame(
        unclassified, columns=["sample", "gene", "zygosity", "reason"])


def normal_filter(df, path):
    ref = read(path)
    if not {'gene', 'category'}.issubset(ref.columns):
        raise ValueError(f'{path} needs gene and category columns')
    if ref['gene'].duplicated().any():
        raise ValueError(f'Duplicate genes in normal reference: {path}')
    ref = ref.set_index('gene')['category']
    # Both conditions required: present in normal reference AND not methylated_on_one.
    present = df['gene'].isin(ref.index)
    categories = df['gene'].map(ref)
    print(f'Absent from normal classification: {(~present).sum()}')
    print(f'Classified methylated_on_one: {(present & categories.eq("methylated_on_one")).sum()}')
    result = df.loc[present & categories.ne('methylated_on_one')].copy()
    result['in_normal'] = result['gene'].map(ref)
    return result


def hgnc_filter(df, path, sample_type):
    hgnc = read(path, low_memory=False)
    if not {'symbol', 'locus_group'}.issubset(hgnc.columns):
        raise ValueError(f'{path} needs symbol and locus_group columns')
    if hgnc['symbol'].duplicated().any():
        raise ValueError(f'Duplicate HGNC symbols in {path}')
    locus = hgnc.set_index('symbol')['locus_group']
    mapping = dict(BASE_MAP)
    pseudo = {'C9orf64'}
    if sample_type == 'BLOOD':
        mapping['THEG'] = 'SPMAP2'
        pseudo.add('C17orf98')
    symbols, groups = [], []
    for gene in df['gene']:
        if gene in locus.index:
            symbol, group = gene, locus.loc[gene]
        elif gene in pseudo:
            symbol, group = gene, 'pseudogene'
        else:
            if gene not in mapping:
                raise ValueError(f'No curated HGNC mapping for {gene!r}')
            symbol = mapping[gene]
            if symbol not in locus.index:
                raise ValueError(f'Mapped symbol {symbol!r} absent from HGNC reference')
            group = locus.loc[symbol]
        symbols.append(symbol)
        groups.append(group)
    result = df.copy()
    result['production_hgnc_gene'] = symbols
    result['locus_group'] = groups
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample-type', choices=['POG', 'BLOOD'], required=True)
    p.add_argument('--data-dir', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True,
                   help='Same output directory used for filter_promoter_DM.py')
    args = p.parse_args()
    root, outdir = args.data_dir, args.output_dir
    paths = sorted((outdir / 'per_sample').glob(f'*{SUFFIX}'))
    if not paths:
        p.error('No per-sample filtered files found; run filter_promoter_DM.py first')
    excluded = set()
    if args.sample_type == 'POG':
        # Reproduce the published analysis: do not recalculate a threshold.
        excluded = set(HYPERMETHYLATED_POG_SAMPLES)
        file_samples = {path.name[:-len(SUFFIX)] for path in paths}
        missing = excluded - file_samples
        if missing:
            print('WARNING: Original excluded samples absent from this run: '
                  + ', '.join(sorted(missing)))
        print(f'Using fixed original POG exclusion list: '
              f'{len(excluded & file_samples)} samples excluded')

        # Audit counts if the filtering stage generated a full-cohort table.
        counts_path = outdir / 'counts/summary_counts_per_sample.tsv'
        if counts_path.is_file():
            counts = read(counts_path)
            required = {'tumour_original_source', 'final_filtered_count'}
            if not required.issubset(counts.columns):
                print('WARNING: Count table lacks columns needed for audit')
            else:
                observed = counts.set_index('tumour_original_source')['final_filtered_count']
                mismatches = {
                    sample: (expected, observed.get(sample, None))
                    for sample, expected in ORIGINAL_HYPERMETHYLATED_COUNTS.items()
                    if sample in observed.index and observed[sample] != expected
                }
                if mismatches:
                    print('WARNING: Original post-imprinted counts differ:')
                    for sample, (expected, actual) in sorted(mismatches.items()):
                        print(f'  {sample}: expected {expected}, observed {actual}')
                else:
                    print('Original hypermethylated-sample counts match '
                          'for all available reference samples')

        exclusion_file = outdir / 'counts/hypermethylated_samples.tsv'
        exclusion_file.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame([
            {'tumour_original_source': sample,
             'original_final_filtered_count': ORIGINAL_HYPERMETHYLATED_COUNTS[sample],
             'present_in_run': sample in file_samples}
            for sample in sorted(excluded)
        ]).to_csv(exclusion_file, sep='\t', index=False)
        print(f'Saved original exclusion list: {exclusion_file}')
    summary, unclassified = summarize(paths, args.sample_type, excluded)
    if args.sample_type == 'POG':
        audit_path = outdir / 'counts/missing_zygosity_observations.tsv'
        audit_path.parent.mkdir(parents=True, exist_ok=True)
        unclassified.to_csv(audit_path, sep='\t', index=False)
        print(f'Unclassified POG DMR observations: {len(unclassified)}; audit: {audit_path}')
    gene_dir = outdir / 'gene_summary'
    save(summary, gene_dir / 'summary_variants_filtered.tsv')
    normal_path = (root / 'annotations/normal_methylation' / args.sample_type /
                   'classification_in_normal_for_genes_w_DMR.tsv.gz')
    df = normal_filter(summary, normal_path)
    save(df, gene_dir / 'summary_after_normal_filter.tsv')
    if args.sample_type == 'POG':
        hpa = read(root / 'annotations/lymphoid_genes/tissue_category_rna_lymphoid_elevated_HPA.tsv.gz')
        if 'Gene' not in hpa:
            raise ValueError('HPA lymphoid reference needs Gene column')
        df = df.loc[~df['gene'].isin(set(hpa['Gene'].dropna()))].copy()
        save(df, gene_dir / 'summary_after_normal_and_lymphoid_filter.tsv')
    annotated = hgnc_filter(df, root / 'annotations/hgnc/hgnc_complete_set.production.tsv.gz', args.sample_type)
    save(annotated, gene_dir / 'summary_with_hgnc_locus_group.tsv')
    protein = annotated.loc[annotated['locus_group'].eq('protein-coding gene')].copy()
    save(protein, gene_dir / 'summary_protein_coding.tsv')
    print(f'Protein-coding genes with total_count > 1: {(protein.total_count > 1).sum()}')
    print(f'Protein-coding genes with total_count > 10: {(protein.total_count > 10).sum()}')


if __name__ == '__main__':
    main()
