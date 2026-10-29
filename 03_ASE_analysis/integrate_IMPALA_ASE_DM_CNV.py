#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Integrate original IMPALA ASE/BAE classifications with promoter DM and CNV.

Python 3.7 compatible. No additional IMPALA adjusted-p-value filtering.
Run from 03_ASE_analysis or specify paths explicitly; see --help.
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

VALID_ASE = ('ASE', 'BAE')
VALID_CNV = ('balance', 'imbalance', 'LOH')
IMPALA_COLUMNS = (
    'gene', 'expression', 'allele1IsMajor', 'majorAlleleFrequency',
    'padj', 'aseResults', 'cnv.A', 'cnv.B', 'cnv_state', 'expectedMAF',
    'methyl_state', 'normalMAF', 'tumorContent'
)
DM_COLUMNS = ('gene', 'pogs_w_aDMs', 'pogs_w_hAMs')
OUTPUT_COLUMNS = (
    'tumour_original_source', 'gene', 'aseResults', 'padj', 'expression',
    'majorAlleleFrequency', 'allele1IsMajor', 'cnv.A', 'cnv.B',
    'cnv_state', 'expectedMAF', 'methyl_state', 'normalMAF',
    'tumorContent', 'dm_status', 'has_aDM', 'has_hDM', 'has_DM',
    'cnv_available', 'cnv_balanced'
)


def parse_args():
    here = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ase-dir',
                        default='/gsc/www/bcgsc.ca/downloads/promoter_methylation_nanopore/ase/IMPALA',
                        help='Published cohort IMPALA directory containing <sample>/summaryTable.tsv.gz')
    parser.add_argument('--dm-file', default=os.path.join(here, 'data', 'summary_protein_coding.tsv'),
                        help='Validated 02_promoter_DM protein-coding summary')
    parser.add_argument('--excluded-samples',
                        default=os.path.join(here, 'data', 'hypermethylated_samples.tsv'),
                        help='DM pipeline hypermethylated-sample TSV with pogid column')
    parser.add_argument('--output-dir', default=os.path.join(here, 'results'),
                        help='Output directory (default: 03_ASE_analysis/results)')
    return parser.parse_args()


def read_tsv(path, required, description):
    if not os.path.isfile(path):
        raise IOError('{} not found: {}'.format(description, path))
    df = pd.read_csv(path, sep='\t', low_memory=False)
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError('{} missing columns: {}'.format(description, ', '.join(missing)))
    return df


def parse_samples(value):
    if pd.isna(value):
        return set()
    return set(s.strip() for s in str(value).split(',') if s.strip())


def build_dm_maps(df):
    adm, hdm = {}, {}
    for _, row in df.iterrows():
        if pd.isna(row['gene']):
            continue
        gene = str(row['gene']).strip()
        if not gene:
            continue
        a = parse_samples(row['pogs_w_aDMs'])
        h = parse_samples(row['pogs_w_hAMs'])
        if a:
            adm.setdefault(gene, set()).update(a)
        if h:
            hdm.setdefault(gene, set()).update(h)
    return adm, hdm


def summarize_sample(df):
    columns = ['aDM', 'hDM', 'no_DM', 'cnv_balance', 'cnv_imbalance', 'cnv_loh']
    result = pd.DataFrame(0, index=list(VALID_ASE), columns=columns, dtype=int)
    for ase in VALID_ASE:
        sub = df.loc[df['aseResults'] == ase]
        result.loc[ase, 'aDM'] = int(sub['has_aDM'].sum())
        result.loc[ase, 'hDM'] = int(sub['has_hDM'].sum())
        result.loc[ase, 'no_DM'] = int((sub['dm_status'] == 'no_DM').sum())
        result.loc[ase, 'cnv_balance'] = int((sub['cnv_state'] == 'balance').sum())
        result.loc[ase, 'cnv_imbalance'] = int((sub['cnv_state'] == 'imbalance').sum())
        result.loc[ase, 'cnv_loh'] = int((sub['cnv_state'] == 'LOH').sum())
    return result


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    count_dir = os.path.join(args.output_dir, 'ASE_BAE_cnvState_count')
    os.makedirs(count_dir, exist_ok=True)

    excluded_df = read_tsv(args.excluded_samples, ['pogid'], 'Exclusion file')
    excluded = set(excluded_df['pogid'].dropna().astype(str).str.strip())
    excluded.discard('')
    if len(excluded) != 14:
        raise ValueError('Expected 14 hypermethylated samples, found {}'.format(len(excluded)))

    dm_df = read_tsv(args.dm_file, DM_COLUMNS, 'Promoter DM file')
    adm_map, hdm_map = build_dm_maps(dm_df)
    print('Cohort defined by deposited IMPALA sample directories')
    print('Hypermethylated exclusions:', len(excluded))
    print('Genes with aDM:', len(adm_map))
    print('Genes with hDM:', len(hdm_map))

    # The published IMPALA directory contains cohort-selected samples.
    # Prefer gzip files if both formats are present for a sample.
    ase_paths = {}
    for pattern in ('summaryTable.tsv', 'summaryTable.tsv.gz'):
        for path in sorted(glob.glob(os.path.join(args.ase_dir, '*', pattern))):
            sample = os.path.basename(os.path.dirname(path))
            if sample in ase_paths and pattern == 'summaryTable.tsv.gz':
                print('INFO: using compressed file for {}'.format(sample))
            ase_paths[sample] = path
    ase_files = [ase_paths[s] for s in sorted(ase_paths)]
    print('IMPALA sample files found:', len(ase_files))
    if not ase_files:
        raise RuntimeError('No IMPALA summaryTable.tsv[.gz] files found under {}'.format(args.ase_dir))
    print('Excluded IMPALA samples present:', len(set(ase_paths).intersection(excluded)))

    gene_tables, sample_summaries, duplicate_tables, processing_stats = [], [], [], []
    seen_samples = set()
    for filename in ase_files:
        sample = os.path.basename(os.path.dirname(filename))
        if sample in excluded:
            continue
        if sample in seen_samples:
            raise RuntimeError('Multiple IMPALA files for sample {}'.format(sample))
        seen_samples.add(sample)
        df = pd.read_csv(filename, sep='\t', low_memory=False)
        missing = [c for c in IMPALA_COLUMNS if c not in df.columns]
        if missing:
            print('WARNING: {} missing {}'.format(sample, ', '.join(missing)), file=sys.stderr)
            processing_stats.append({'sample': sample, 'status': 'missing_columns',
                                     'missing_columns': ','.join(missing),
                                     'original_genes': len(df)})
            continue
        original_n = len(df)
        df = df.loc[df['aseResults'].isin(VALID_ASE) & df['gene'].notna()].copy()
        df['gene'] = df['gene'].astype(str).str.strip()
        df = df.loc[df['gene'] != ''].copy()

        duplicate_mask = df.duplicated(subset=['gene'], keep=False)
        duplicates = df.loc[duplicate_mask].copy()
        if not duplicates.empty:
            duplicates.insert(0, 'tumour_original_source', sample)
            duplicate_tables.append(duplicates)
            df = df.loc[~duplicate_mask].copy()

        df['tumour_original_source'] = sample
        df['has_aDM'] = [sample in adm_map.get(g, set()) for g in df['gene']]
        df['has_hDM'] = [sample in hdm_map.get(g, set()) for g in df['gene']]
        df['has_DM'] = df['has_aDM'] | df['has_hDM']
        df['dm_status'] = np.select(
            [df['has_aDM'] & df['has_hDM'], df['has_aDM'], df['has_hDM']],
            ['both_aDM_hDM', 'aDM', 'hDM'], default='no_DM')
        df['cnv_state'] = df['cnv_state'].where(df['cnv_state'].isin(VALID_CNV), np.nan)
        df['cnv_available'] = df['cnv_state'].notna()
        df['cnv_balanced'] = df['cnv_state'] == 'balance'
        output = df.loc[:, list(OUTPUT_COLUMNS)].copy()
        gene_tables.append(output)

        summary = summarize_sample(output)
        summary.to_csv(os.path.join(count_dir, '{}_ASE_BAE_cnvState_count.tsv'.format(sample)), sep='\t')
        summary_long = summary.reset_index().rename(columns={'index': 'aseResults'})
        summary_long.insert(0, 'tumour_original_source', sample)
        sample_summaries.append(summary_long)
        processing_stats.append({
            'sample': sample, 'status': 'processed', 'missing_columns': '',
            'original_genes': original_n, 'retained_genes': len(output),
            'duplicate_rows_excluded': len(duplicates),
            'ASE': int((output['aseResults'] == 'ASE').sum()),
            'BAE': int((output['aseResults'] == 'BAE').sum()),
            'aDM': int(output['has_aDM'].sum()),
            'hDM': int(output['has_hDM'].sum()),
            'CNV_available': int(output['cnv_available'].sum()),
            'methyl_available': int(output['methyl_state'].notna().sum())
        })

    if not gene_tables:
        raise RuntimeError('No usable IMPALA files remained after exclusions')
    combined = pd.concat(gene_tables, ignore_index=True)
    combined.to_csv(os.path.join(args.output_dir, 'IMPALA_ASE_DM_CNV_gene_sample_table.tsv'),
                    sep='\t', index=False)
    pd.concat(sample_summaries, ignore_index=True).to_csv(
        os.path.join(args.output_dir, 'IMPALA_ASE_BAE_sample_summary.tsv'), sep='\t', index=False)
    if duplicate_tables:
        pd.concat(duplicate_tables, ignore_index=True).to_csv(
            os.path.join(args.output_dir, 'IMPALA_duplicate_gene_records.tsv'), sep='\t', index=False)
    else:
        pd.DataFrame(columns=['tumour_original_source'] + list(IMPALA_COLUMNS)).to_csv(
            os.path.join(args.output_dir, 'IMPALA_duplicate_gene_records.tsv'), sep='\t', index=False)
    pd.DataFrame(processing_stats).to_csv(
        os.path.join(args.output_dir, 'IMPALA_processing_summary.tsv'), sep='\t', index=False)

    print('\n' + '=' * 70 + '\nFINAL SUMMARY\n' + '=' * 70)
    print('Samples processed:', combined['tumour_original_source'].nunique())
    print('Gene-sample observations:', len(combined))
    print('\nASE classifications:\n{}'.format(combined['aseResults'].value_counts()))
    print('\nDM classifications:\n{}'.format(combined['dm_status'].value_counts()))
    print('\nCopy-number classifications:\n{}'.format(combined['cnv_state'].value_counts(dropna=False)))
    print('\nASE versus DM:\n{}'.format(pd.crosstab(combined['aseResults'], combined['dm_status'])))
    print('\nASE versus copy-number state:\n{}'.format(
        pd.crosstab(combined['aseResults'], combined['cnv_state'], dropna=False)))
    balanced = combined.loc[combined['cnv_state'] == 'balance']
    print('\nASE and aDM in balanced CN regions:\n{}'.format(
        pd.crosstab(balanced['aseResults'], balanced['has_aDM'])))
    print('\nOutput directory:', args.output_dir)
    print('Analysis complete.')


if __name__ == '__main__':
    main()
