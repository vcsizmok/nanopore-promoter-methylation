#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Directional concordance between promoter aDM and ASE (Python 3.7 compatible).

Input: integrated, cohort-filtered IMPALA gene-sample table produced by
integrate_IMPALA_ASE_DM_CNV.py. No further sample or padj filtering.

Concordance: methyl_state = HP1 - HP2; allele1IsMajor=True means HP1
is more highly expressed. Concordant when the less methylated haplotype
is the more highly expressed haplotype.

The gene-sample binomial test assumes independent observations; use the
sample-level sign test as a complementary descriptive analysis.
"""
from __future__ import print_function

import argparse
import os
import numpy as np
import pandas as pd
from scipy.stats import beta
try:
    from scipy.stats import binomtest
    def binomial_p(k, n):
        return float(binomtest(k, n, 0.5, alternative='greater').pvalue)
except ImportError:
    from scipy.stats import binom_test
    def binomial_p(k, n):
        return float(binom_test(k, n, 0.5, alternative='greater'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_INPUT = os.path.join(SCRIPT_DIR, 'results', 'IMPALA_ASE_DM_CNV_gene_sample_table.tsv')
DEFAULT_OUTPUT = os.path.join(SCRIPT_DIR, 'results', 'directional_aDM_ASE_concordance')


def bool_parse(value):
    if pd.isnull(value):
        return np.nan
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    value = str(value).strip().lower()
    if value in ('true', 't', '1', '1.0'):
        return True
    if value in ('false', 'f', '0', '0.0'):
        return False
    return np.nan


def exact_ci(k, n, alpha=0.05):
    if n == 0:
        return np.nan, np.nan
    low = 0.0 if k == 0 else beta.ppf(alpha / 2.0, k, n - k + 1)
    high = 1.0 if k == n else beta.ppf(1.0 - alpha / 2.0, k + 1, n - k)
    return float(low), float(high)


def summarize(df, label):
    n = len(df)
    k = int(df['less_methylated_more_expressed'].sum()) if n else 0
    low, high = exact_ci(k, n)
    return dict(group=label, n_informative=n, n_concordant=k,
                n_discordant=n-k,
                concordance_percent=100.0*k/n if n else np.nan,
                exact_CI_low_percent=100.0*low,
                exact_CI_high_percent=100.0*high,
                binomial_one_sided_p=binomial_p(k, n) if n else np.nan,
                n_samples=df['tumour_original_source'].nunique(),
                n_genes=df['gene'].nunique())


def sample_summary(df, label):
    columns = ['tumour_original_source', 'n_events', 'n_concordant',
               'n_discordant', 'concordance_percent', 'group']
    if df.empty:
        return pd.DataFrame(columns=columns)
    result = (df.groupby('tumour_original_source')['less_methylated_more_expressed']
              .agg(['size', 'sum']).reset_index())
    result.rename(columns={'size': 'n_events', 'sum': 'n_concordant'}, inplace=True)
    result['n_concordant'] = result['n_concordant'].astype(int)
    result['n_discordant'] = result['n_events'] - result['n_concordant']
    result['concordance_percent'] = 100.0 * result['n_concordant'] / result['n_events']
    result['group'] = label
    return result[columns]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', default=DEFAULT_INPUT,
                        help='Integrated IMPALA gene-sample TSV (default: results/...)')
    parser.add_argument('--output-dir', default=DEFAULT_OUTPUT)
    parser.add_argument('--min-abs-methyl', type=float, default=0.20)
    parser.add_argument('--min-sample-events', type=int, default=5)
    args = parser.parse_args()
    if args.min_abs_methyl < 0 or args.min_sample_events < 1:
        parser.error('Thresholds must be nonnegative and min-sample-events >= 1')
    if not os.path.isfile(args.input):
        raise IOError('Input not found: ' + args.input)
    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)

    df = pd.read_csv(args.input, sep='\t', low_memory=False)
    required = ['tumour_original_source', 'gene', 'aseResults', 'dm_status',
                'cnv_state', 'methyl_state', 'allele1IsMajor']
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError('Missing columns: ' + ', '.join(missing))

    d = df.loc[(df['aseResults'] == 'ASE') & (df['dm_status'] == 'aDM')].copy()
    d['methyl_difference_HP1_minus_HP2'] = pd.to_numeric(d['methyl_state'], errors='coerce')
    d['HP1_more_expressed'] = d['allele1IsMajor'].map(bool_parse)
    d['has_methyl_difference'] = d['methyl_difference_HP1_minus_HP2'].notnull()
    d['has_expression_direction'] = d['HP1_more_expressed'].notnull()
    d['nonzero_methyl_difference'] = (d['has_methyl_difference'] &
                                       (d['methyl_difference_HP1_minus_HP2'] != 0))
    d['direction_informative'] = d['nonzero_methyl_difference'] & d['has_expression_direction']

    # Object dtype avoids nullable-boolean incompatibilities in older pandas.
    d['less_methylated_more_expressed'] = pd.Series([np.nan] * len(d), index=d.index, dtype=object)
    good = d['direction_informative']
    if good.any():
        hp1_less_methylated = (d.loc[good, 'methyl_difference_HP1_minus_HP2'] < 0).values
        hp1_more_expressed = d.loc[good, 'HP1_more_expressed'].astype(bool).values
        d.loc[good, 'less_methylated_more_expressed'] = (
            hp1_less_methylated == hp1_more_expressed).tolist()
    d['abs_methyl_difference'] = d['methyl_difference_HP1_minus_HP2'].abs()
    d.to_csv(os.path.join(args.output_dir, 'aDM_ASE_directional_all_observations.tsv'),
             sep='\t', index=False)

    audit = []
    for state in ('all', 'balance', 'imbalance', 'LOH'):
        sub = d if state == 'all' else d.loc[d['cnv_state'] == state]
        audit.append(dict(cnv_group=state, n_aDM_ASE=len(sub),
                          missing_methyl=int((~sub['has_methyl_difference']).sum()),
                          missing_expression_direction=int((~sub['has_expression_direction']).sum()),
                          zero_methyl=int((sub['has_methyl_difference'] &
                                           ~sub['nonzero_methyl_difference']).sum()),
                          n_direction_informative=int(sub['direction_informative'].sum())))
    audit_df = pd.DataFrame(audit)
    audit_df.to_csv(os.path.join(args.output_dir, 'directional_data_availability_by_CNV.tsv'),
                    sep='\t', index=False)

    informative = d.loc[d['direction_informative']].copy()
    informative['less_methylated_more_expressed'] = informative['less_methylated_more_expressed'].astype(bool)
    informative.to_csv(os.path.join(args.output_dir, 'aDM_ASE_directional_informative.tsv'),
                       sep='\t', index=False)

    balanced = informative.loc[informative['cnv_state'] == 'balance'].copy()
    balanced_sample = sample_summary(balanced, 'balanced_primary')
    eligible = set(balanced_sample.loc[balanced_sample['n_events'] >= args.min_sample_events,
                                       'tumour_original_source'])
    groups = [
        ('balanced_primary', balanced),
        ('all_CNV_secondary', informative),
        ('imbalanced_secondary', informative.loc[informative['cnv_state'] == 'imbalance']),
        ('LOH_descriptive', informative.loc[informative['cnv_state'] == 'LOH']),
        ('balanced_abs_methyl_ge_{:g}'.format(args.min_abs_methyl),
         balanced.loc[balanced['abs_methyl_difference'] >= args.min_abs_methyl]),
        ('balanced_samples_ge{}_informative_events'.format(args.min_sample_events),
         balanced.loc[balanced['tumour_original_source'].isin(eligible)])
    ]
    per_sample = pd.concat([sample_summary(sub, label) for label, sub in groups],
                           ignore_index=True)
    per_sample.to_csv(os.path.join(args.output_dir, 'per_sample_directional_concordance.tsv'),
                      sep='\t', index=False)
    summary_df = pd.DataFrame([summarize(sub, label) for label, sub in groups])
    summary_df.to_csv(os.path.join(args.output_dir, 'directional_concordance_summary.tsv'),
                      sep='\t', index=False)

    signs = []
    for label, sub in groups:
        if label not in ('balanced_primary',
                         'balanced_samples_ge{}_informative_events'.format(args.min_sample_events)):
            continue
        ss = sample_summary(sub, label)
        positive = int((ss['concordance_percent'] > 50).sum())
        negative = int((ss['concordance_percent'] < 50).sum())
        ties = int((ss['concordance_percent'] == 50).sum())
        signs.append(dict(group=label, n_samples=len(ss), samples_above_50=positive,
                          samples_below_50=negative, samples_equal_50=ties,
                          sign_test_one_sided_p=(binomial_p(positive, positive+negative)
                                                 if positive+negative else np.nan)))
    sign_df = pd.DataFrame(signs)
    sign_df.to_csv(os.path.join(args.output_dir, 'sample_level_sign_tests.tsv'),
                   sep='\t', index=False)

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    if not balanced_sample.empty:
        ax.scatter(balanced_sample['n_events'], balanced_sample['concordance_percent'],
                   alpha=0.55, s=16)
        ax.set_xscale('log')
    ax.axhline(50, linestyle='--', linewidth=1, color='gray')
    ax.set_xlabel('Informative aDM + ASE observations per sample (log scale)')
    ax.set_ylabel('Less-methylated allele more expressed (%)')
    ax.set_ylim(-3, 103)
    fig.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(args.output_dir,
                                 'balanced_per_sample_concordance_vs_event_count.' + ext), dpi=300)
    plt.close(fig)

    print('\nDATA AVAILABILITY\n' + audit_df.to_string(index=False))
    print('\nDIRECTIONAL CONCORDANCE\n' + summary_df.to_string(index=False))
    print('\nSAMPLE-LEVEL SIGN TESTS\n' + sign_df.to_string(index=False))
    print('\nOutput: ' + args.output_dir)
    print('NOTE: Gene-sample binomial tests do not account for within-sample or recurrent-gene dependence.')


if __name__ == '__main__':
    main()
