#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compare ASE versus BAE frequencies of promoter aDM and CN states.

Input is the cohort-filtered gene-sample table produced by
integrate_IMPALA_ASE_DM_CNV.py. No additional cohort exclusions or
IMPALA adjusted-p-value filters are applied.

Python 3.7 compatible.
"""
from __future__ import print_function

import argparse
import os
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, fisher_exact

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

METRICS = [
    ('aDM', 'ASE_aDM', 'BAE_aDM'),
    ('LOH', 'ASE_LOH', 'BAE_LOH'),
    ('Imbalanced CN', 'ASE_Imbalance', 'BAE_Imbalance'),
    ('Balanced CN', 'ASE_Balance', 'BAE_Balance'),
]


def paired_wilcoxon(x, y):
    valid = np.isfinite(x) & np.isfinite(y)
    a, b = x[valid], y[valid]
    if len(a) == 0:
        return np.nan, 0
    if np.all(a == b):
        return 1.0, len(a)
    return float(wilcoxon(a, b, zero_method='zsplit').pvalue), len(a)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', default=os.path.join(here, 'results',
                        'IMPALA_ASE_DM_CNV_gene_sample_table.tsv'))
    parser.add_argument('--output-dir', default=os.path.join(here, 'results',
                        'ASE_BAE_cnvState_count', 'wilcoxon_unfiltered'))
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    required = ['tumour_original_source', 'gene', 'aseResults', 'dm_status', 'cnv_state']
    df = pd.read_csv(args.input, sep='\t', usecols=required, low_memory=False)
    for col in required:
        df[col] = df[col].astype(str).str.strip()
    df = df[df['aseResults'].isin(['ASE', 'BAE'])].copy()
    if df.empty:
        raise ValueError('No ASE/BAE observations in input')
    if df.duplicated(['tumour_original_source', 'gene']).any():
        raise ValueError('Duplicate gene-sample records in integrated input')

    df['is_aDM'] = df['dm_status'].isin(['aDM', 'both_aDM_hDM'])
    df['is_hDM'] = df['dm_status'].eq('hDM')
    df['is_noDM'] = df['dm_status'].eq('no_DM')
    df['valid_dm'] = df['is_aDM'] | df['is_hDM'] | df['is_noDM']
    df['valid_cn'] = df['cnv_state'].isin(['balance', 'imbalance', 'LOH'])

    print('Samples in integrated input:', df['tumour_original_source'].nunique())
    print('ASE/BAE observations:', len(df))
    print('Unknown DM classifications:', int((~df['valid_dm']).sum()))
    print('Missing/unrecognized CN states:', int((~df['valid_cn']).sum()))

    for source, target in [('is_aDM', 'aDM_count'), ('is_hDM', 'hDM_count'),
                           ('is_noDM', 'noDM_count'), ('valid_dm', 'dm_denom'),
                           ('valid_cn', 'cn_denom')]:
        df[target] = df[source].astype(int)
    for state, target in [('balance', 'cnv_balance'), ('imbalance', 'cnv_imbalance'),
                          ('LOH', 'cnv_loh')]:
        df[target] = df['cnv_state'].eq(state).astype(int)
    count_cols = ['aDM_count', 'hDM_count', 'noDM_count', 'cnv_balance',
                  'cnv_imbalance', 'cnv_loh', 'dm_denom', 'cn_denom']
    grouped = df.groupby(['tumour_original_source', 'aseResults'])[count_cols].sum()
    grouped.to_csv(os.path.join(args.output_dir, 'ASE_BAE_per_sample_counts.tsv'), sep='\t')
    totals = df.groupby('aseResults')[count_cols].sum()
    totals.to_csv(os.path.join(args.output_dir, 'ASE_BAE_pooled_counts.tsv'), sep='\t')

    records = []
    for sample in sorted(df['tumour_original_source'].unique()):
        rec = {'sample': sample}
        for label in ['ASE', 'BAE']:
            key = (sample, label)
            # An absent ASE/BAE class must not be mistaken for zero percent.
            if key not in grouped.index:
                for name in ['n_DM_known', 'n_CNV_known', 'aDM', 'LOH',
                             'Imbalance', 'Balance']:
                    rec[label + '_' + name] = np.nan
                continue
            row = grouped.loc[key]
            dm_total, cn_total = float(row['dm_denom']), float(row['cn_denom'])
            rec[label + '_n_DM_known'] = int(dm_total)
            rec[label + '_n_CNV_known'] = int(cn_total)
            rec[label + '_aDM'] = 100.0 * row['aDM_count'] / dm_total if dm_total else np.nan
            for name, col in [('LOH', 'cnv_loh'), ('Imbalance', 'cnv_imbalance'),
                              ('Balance', 'cnv_balance')]:
                rec[label + '_' + name] = 100.0 * row[col] / cn_total if cn_total else np.nan
        records.append(rec)
    percentages = pd.DataFrame(records)
    percentages.to_csv(os.path.join(args.output_dir, 'ASE_BAE_aDM_cnv_percentages.tsv'),
                       sep='\t', index=False)

    stats_rows = []
    for title, ase_col, bae_col in METRICS:
        x = percentages[ase_col].values.astype(float)
        y = percentages[bae_col].values.astype(float)
        p, n = paired_wilcoxon(x, y)
        valid = np.isfinite(x) & np.isfinite(y)
        diff = x[valid] - y[valid]
        stats_rows.append({
            'comparison': title, 'n_paired_samples': n,
            'median_ASE_percent': float(np.median(x[valid])) if n else np.nan,
            'median_BAE_percent': float(np.median(y[valid])) if n else np.nan,
            'median_paired_difference_pp': float(np.median(diff)) if n else np.nan,
            'samples_ASE_greater': int(np.sum(diff > 0)),
            'samples_BAE_greater': int(np.sum(diff < 0)),
            'samples_tied': int(np.sum(diff == 0)),
            'wilcoxon_signed_rank_zsplit_p': p,
        })
    stats_df = pd.DataFrame(stats_rows)
    stats_df.to_csv(os.path.join(args.output_dir, 'ASE_BAE_paired_Wilcoxon_results.tsv'),
                    sep='\t', index=False)
    print('\nPAIRED WILCOXON RESULTS\n' + stats_df.to_string(index=False))

    if 'ASE' in totals.index and 'BAE' in totals.index:
        tab = [[int(totals.loc['ASE', 'aDM_count']),
                int(totals.loc['ASE', 'hDM_count'] + totals.loc['ASE', 'noDM_count'])],
               [int(totals.loc['BAE', 'aDM_count']),
                int(totals.loc['BAE', 'hDM_count'] + totals.loc['BAE', 'noDM_count'])]]
        odds, p = fisher_exact(tab)
        pd.DataFrame([{'ASE_aDM': tab[0][0], 'ASE_other': tab[0][1],
                       'BAE_aDM': tab[1][0], 'BAE_other': tab[1][1],
                       'pooled_fisher_OR': odds, 'pooled_fisher_p': p}]).to_csv(
            os.path.join(args.output_dir, 'aDM_pooled_Fisher_descriptive.tsv'),
            sep='\t', index=False)

    fig, axes = plt.subplots(1, 4, figsize=(14, 4), dpi=150)
    for ax, (title, ase_col, bae_col), stat in zip(axes, METRICS, stats_rows):
        x = percentages[ase_col].dropna().values
        y = percentages[bae_col].dropna().values
        values = [x, y]
        if len(x) and len(y):
            violin = ax.violinplot(values, positions=[1, 2], showextrema=False)
            for body in violin['bodies']:
                body.set_alpha(0.35)
            ax.boxplot(values, positions=[1, 2], widths=0.3, showfliers=False,
                       medianprops={'color': 'black', 'linewidth': 1.5})
        ax.set_xticks([1, 2])
        ax.set_xticklabels(['ASE', 'BAE'])
        ax.set_title('{}\npaired n={}, p={:.2e}'.format(
            title, stat['n_paired_samples'], stat['wilcoxon_signed_rank_zsplit_p']), fontsize=10)
        ax.set_ylabel('% of genes' if ax is axes[0] else '')
        ax.set_ylim(0, 100)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    fig.tight_layout()
    for ext in ['png', 'pdf']:
        fig.savefig(os.path.join(args.output_dir, 'ASE_BAE_aDM_CNV_unfiltered.' + ext),
                    dpi=300, bbox_inches='tight')
    plt.close(fig)
    print('\nOutputs:', args.output_dir)
    print('NOTE: Pooled Fisher test is descriptive and does not account for sample clustering.')


if __name__ == '__main__':
    main()
