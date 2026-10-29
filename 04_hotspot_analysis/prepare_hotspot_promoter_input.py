#!/usr/bin/env python3
"""Prepare the promoter universe and binary DM labels for hotspot analysis.

Uses the refTSS single-best-transcript annotation and the validated
02_promoter_DM protein-coding summary shared with 03_ASE_analysis.
Python 3.7 compatible.
"""
import argparse
import os
import pandas as pd
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_TSS = '/gsc/www/bcgsc.ca/downloads/promoter_methylation_nanopore/annotations/refTSS/refTSS_hg38_single_best_transcript.bed.gz'
DEFAULT_DM = os.path.normpath(os.path.join(HERE, '..', '03_ASE_analysis', 'data', 'summary_protein_coding.tsv'))
DEFAULT_OUT = os.path.join(HERE, 'data', 'input_promoters_coordinates_methylation.tsv')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reftss', default=DEFAULT_TSS, help='refTSS BED-like headered TSV')
    parser.add_argument('--dm-summary', default=DEFAULT_DM, help='Validated protein-coding DM summary')
    parser.add_argument('--output', default=DEFAULT_OUT)
    parser.add_argument('--upstream', type=int, default=1500)
    parser.add_argument('--downstream', type=int, default=1000)
    args = parser.parse_args()
    if args.upstream < 0 or args.downstream < 0:
        parser.error('upstream/downstream must be nonnegative')

    tss = pd.read_csv(args.reftss, sep='\t', low_memory=False)
    dm = pd.read_csv(args.dm_summary, sep='\t', low_memory=False)
    required = {'chrom', 'start', 'strand', 'symbol'}
    if not required.issubset(tss.columns):
        raise ValueError('refTSS missing columns: ' + ', '.join(sorted(required - set(tss.columns))))
    if 'gene' not in dm.columns:
        raise ValueError('DM summary missing gene column')
    if tss['symbol'].isna().any() or dm['gene'].isna().any():
        raise ValueError('Missing gene identifiers in input')
    if tss['symbol'].duplicated().any():
        examples = tss.loc[tss['symbol'].duplicated(False), 'symbol'].unique()[:10]
        raise ValueError('Duplicate refTSS symbols (original script assumes unique symbols): ' + ', '.join(map(str, examples)))
    if dm['gene'].duplicated().any():
        raise ValueError('Duplicate gene identifiers in DM summary')

    tss = tss.loc[~tss['chrom'].isin(['chrX', 'chrY', 'chrM'])].copy()
    if not tss['strand'].isin(['+', '-']).all():
        raise ValueError('Unexpected strand values')
    tss['start'] = pd.to_numeric(tss['start'], errors='raise')
    if tss['start'].isna().any() or (tss['start'] % 1 != 0).any():
        raise ValueError('Invalid refTSS start coordinates')
    tss['start'] = tss['start'].astype(np.int64)

    # Original code uses the refTSS start coordinate as TSS on both strands.
    # Normalize interval boundaries for negative-strand genes; midpoint is unchanged.
    plus = tss['strand'].eq('+')
    a = np.where(plus, tss['start'] - args.upstream, tss['start'] + args.upstream)
    b = np.where(plus, tss['start'] + args.downstream, tss['start'] - args.downstream)
    tss['promoter_start'] = np.minimum(a, b).astype(np.int64)
    tss['promoter_end'] = np.maximum(a, b).astype(np.int64)
    if (tss['promoter_start'] < 0).any():
        raise ValueError('Negative promoter start; review chromosome-boundary handling')
    tss['methylated'] = tss['symbol'].isin(set(dm['gene'])).astype(int)
    tss = tss.rename(columns={'start': 'tss_start'})
    output = tss[['symbol', 'chrom', 'tss_start', 'strand', 'promoter_start', 'promoter_end', 'methylated']]
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    output.to_csv(args.output, sep='\t', index=False)
    print('Autosomal refTSS promoters: {:,}'.format(len(output)))
    print('DM-positive promoters: {:,}'.format(int(output['methylated'].sum())))
    print('DM-negative promoters: {:,}'.format(int((output['methylated'] == 0).sum())))
    print('Output: ' + args.output)
    if output['methylated'].nunique() < 2:
        print('WARNING: Only one DM class present; verify DM summary contains only DM-positive genes.')


if __name__ == '__main__':
    main()
