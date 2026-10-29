#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Sample-stratified association of promoter aDM with ASE across CN states.

Primary: balanced copy number, aDM versus no_DM.
Secondary: imbalance, LOH, all known CN states.
Sensitivity: balanced samples with >=5 aDM genes.

Uses the integrated, unfiltered IMPALA ASE/BAE classifications; no additional
padj filtering or hypermethylated-sample exclusions are applied here.
Python 3.7 compatible.
"""
from __future__ import print_function
import argparse
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.contingency_tables import StratifiedTable
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
GROUPS = ('balanced', 'imbalanced', 'LOH', 'all_known_CNV')
CELLS = ('aDM_ASE', 'aDM_BAE', 'noDM_ASE', 'noDM_BAE')


def sign_p(pos, neg):
    if pos + neg == 0:
        return np.nan
    if hasattr(stats, 'binomtest'):
        return float(stats.binomtest(pos, pos + neg, p=0.5, alternative='two-sided').pvalue)
    return float(stats.binom_test(pos, pos + neg, p=0.5, alternative='two-sided'))


def sample_counts(df):
    records = []
    for (sample, cn), sub in df.groupby(['sample', 'cnv_group']):
        counts = dict((c, 0) for c in CELLS)
        for dm, ase in zip(sub['dm_status'], sub['aseResults']):
            counts[('aDM' if dm == 'aDM' else 'noDM') + '_' + ase] += 1
        a, b, c, d = [counts[k] for k in CELLS]
        n_adm, n_nodm = a + b, c + d
        adm_pct = 100.0 * a / n_adm if n_adm else np.nan
        nodm_pct = 100.0 * c / n_nodm if n_nodm else np.nan
        odds = ((a + .5) * (d + .5)) / ((b + .5) * (c + .5))
        # Both row and column margins must be nonzero for a CMH stratum.
        informative = n_adm > 0 and n_nodm > 0 and (a + c) > 0 and (b + d) > 0
        rec = dict(sample=sample, cnv_group=cn, aDM_total=n_adm,
                   noDM_total=n_nodm, ASE_percent_aDM=adm_pct,
                   ASE_percent_noDM=nodm_pct,
                   percentage_point_difference=adm_pct - nodm_pct,
                   odds_ratio_corrected=odds, log2_odds_ratio=np.log2(odds),
                   informative_for_CMH=informative)
        rec.update(counts)
        records.append(rec)
    return pd.DataFrame(records)


def summarize(counts, group, analysis, min_adm=1):
    sub = counts[counts['cnv_group'] == group].copy()
    if analysis == 'balanced_aDM_ge5':
        sub = sub[(sub['aDM_total'] >= min_adm) & sub['informative_for_CMH']].copy()
    informative = sub[sub['informative_for_CMH']].copy()
    rec = dict(analysis=analysis, cnv_group=group,
               min_aDM_genes=min_adm if analysis == 'balanced_aDM_ge5' else 1,
               n_samples=len(sub), n_informative_samples=len(informative))
    for col in CELLS:
        rec['n_' + col] = int(sub[col].sum())
    rec.update(CMH_OR=np.nan, CMH_CI_low=np.nan, CMH_CI_high=np.nan,
               CMH_p=np.nan, heterogeneity_p=np.nan, median_sample_OR=np.nan,
               positive_samples=0, negative_samples=0, tied_samples=0,
               sign_test_p=np.nan)
    if not informative.empty:
        diffs = informative['percentage_point_difference']
        rec['positive_samples'] = int((diffs > 1e-10).sum())
        rec['negative_samples'] = int((diffs < -1e-10).sum())
        rec['tied_samples'] = int((diffs.abs() <= 1e-10).sum())
        rec['sign_test_p'] = sign_p(rec['positive_samples'], rec['negative_samples'])
        rec['median_sample_OR'] = float(informative['odds_ratio_corrected'].median())
    if len(informative) >= 2:
        tables = [np.array([[r.aDM_ASE, r.aDM_BAE],
                            [r.noDM_ASE, r.noDM_BAE]], dtype=float)
                  for r in informative.itertuples(index=False)]
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            st = StratifiedTable(tables)
            rec['CMH_OR'] = float(st.oddsratio_pooled)
            lo, hi = st.oddsratio_pooled_confint()
            rec['CMH_CI_low'], rec['CMH_CI_high'] = float(lo), float(hi)
            rec['CMH_p'] = float(st.test_null_odds().pvalue)
            try:
                rec['heterogeneity_p'] = float(st.test_equal_odds().pvalue)
            except (ValueError, ZeroDivisionError):
                pass
    return rec


def save_fig(fig, outdir, stem):
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(str(outdir / (stem + '.' + ext)), dpi=300, bbox_inches='tight')
    plt.close(fig)


def plot_paired(counts, outdir, stem, title):
    sub = counts[counts['informative_for_CMH']].dropna(
        subset=['ASE_percent_aDM', 'ASE_percent_noDM'])
    if sub.empty:
        return
    fig, ax = plt.subplots(figsize=(5.5, 5))
    for row in sub.itertuples(index=False):
        ax.plot([0, 1], [row.ASE_percent_noDM, row.ASE_percent_aDM],
                color='0.6', alpha=.17, linewidth=.55)
    ax.boxplot([sub['ASE_percent_noDM'].values, sub['ASE_percent_aDM'].values],
               positions=[0, 1], widths=.35, showfliers=False)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['No promoter DM', 'aDM'])
    ax.set_ylabel('Genes exhibiting ASE (%)')
    ax.set_ylim(-3, 103)
    ax.set_title('{} (n = {} samples)'.format(title, len(sub)))
    save_fig(fig, outdir, stem)


def plot_forest(summary, outdir):
    sub = summary[(summary['analysis'] == 'primary') &
                  np.isfinite(summary['CMH_OR']) &
                  np.isfinite(summary['CMH_CI_low']) &
                  np.isfinite(summary['CMH_CI_high']) &
                  (summary['CMH_CI_low'] > 0)].copy()
    if sub.empty:
        return
    fig, ax = plt.subplots(figsize=(6.5, 3.5))
    y = np.arange(len(sub))
    odds = sub['CMH_OR'].values
    lo = sub['CMH_CI_low'].values
    hi = sub['CMH_CI_high'].values
    ax.errorbar(odds, y, xerr=[np.maximum(0, odds-lo), np.maximum(0, hi-odds)],
                fmt='o', capsize=3)
    ax.axvline(1, color='0.5', linestyle='--')
    ax.set_xscale('log')
    ax.set_yticks(y)
    ax.set_yticklabels(sub['cnv_group'])
    ax.set_xlabel('CMH odds ratio (95% CI)')
    ax.invert_yaxis()
    save_fig(fig, outdir, 'CNV_stratified_CMH_forest')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', default=str(ROOT / 'results' /
                        'IMPALA_ASE_DM_CNV_gene_sample_table.tsv'))
    parser.add_argument('--output-dir', default=str(ROOT / 'results' /
                        'copy_number_stratified_aDM_ASE'))
    parser.add_argument('--min-adm', type=int, default=5,
                        help='Minimum balanced-CN aDM genes per sample for sensitivity (default: 5)')
    args = parser.parse_args()
    if args.min_adm < 1:
        parser.error('--min-adm must be >= 1')
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.input, sep='\t', low_memory=False)
    required = ['tumour_original_source', 'gene', 'aseResults', 'dm_status',
                'cnv_state', 'majorAlleleFrequency', 'expectedMAF']
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError('Missing columns: {}'.format(', '.join(missing)))
    df = df[df['dm_status'].isin(['aDM', 'no_DM']) &
            df['aseResults'].isin(['ASE', 'BAE'])].copy()
    df['sample'] = df['tumour_original_source'].astype(str).str.strip()
    df['cnv_group'] = df['cnv_state'].astype(str).str.strip().str.lower().map(
        {'balance': 'balanced', 'balanced': 'balanced',
         'imbalance': 'imbalanced', 'imbalanced': 'imbalanced', 'loh': 'LOH'})
    duplicate = df.duplicated(['sample', 'gene'], keep=False)
    if duplicate.any():
        raise ValueError('{} duplicate gene-sample rows; resolve upstream'.format(int(duplicate.sum())))
    known = df[df['cnv_group'].notna()].copy()
    if known.empty:
        raise ValueError('No rows with recognized copy-number states')
    all_known = known.copy()
    all_known['cnv_group'] = 'all_known_CNV'
    counts = sample_counts(pd.concat([known, all_known], ignore_index=True))
    counts['eligible_balanced_aDM_ge5'] = ((counts['cnv_group'] == 'balanced') &
        (counts['aDM_total'] >= args.min_adm) & counts['informative_for_CMH'])
    counts.to_csv(str(outdir / 'sample_specific_aDM_ASE_results.tsv'), sep='\t', index=False)
    results = [summarize(counts, g, 'primary') for g in GROUPS]
    results.append(summarize(counts, 'balanced', 'balanced_aDM_ge5', args.min_adm))
    summary = pd.DataFrame(results)
    summary.to_csv(str(outdir / 'CMH_analysis_summary_with_sensitivity.tsv'),
                   sep='\t', index=False)
    summary[summary['analysis'] == 'primary'].drop(
        columns=['analysis', 'min_aDM_genes']).to_csv(
            str(outdir / 'CMH_analysis_summary.tsv'), sep='\t', index=False)
    sensitivity = counts[counts['eligible_balanced_aDM_ge5']].copy()
    sensitivity.to_csv(str(outdir / 'balanced_CNV_aDM_ge5_sample_results.tsv'),
                       sep='\t', index=False)
    summary[summary['analysis'] == 'balanced_aDM_ge5'].to_csv(
        str(outdir / 'balanced_CNV_aDM_ge5_CMH_summary.tsv'), sep='\t', index=False)
    balanced = counts[counts['cnv_group'] == 'balanced']
    plot_paired(balanced, outdir, 'balanced_CNV_paired_ASE_frequency',
                'Balanced CN: all informative')
    plot_paired(sensitivity, outdir, 'balanced_CNV_aDM_ge5_paired_ASE_frequency',
                'Balanced CN: >= {} aDM genes'.format(args.min_adm))
    plot_forest(summary, outdir)
    hist = balanced[balanced['informative_for_CMH']]['log2_odds_ratio'].replace(
        [np.inf, -np.inf], np.nan).dropna()
    if not hist.empty:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.hist(hist, bins=35)
        ax.axvline(0, linestyle='--', color='0.4')
        ax.set_xlabel('Sample-specific log2 odds ratio (+0.5 correction)')
        ax.set_ylabel('Number of samples')
        save_fig(fig, outdir, 'balanced_CNV_sample_OR_histogram')
    # Preserve the original observed-minus-expected major allele frequency audit.
    maf = df.copy()
    maf['majorAlleleFrequency'] = pd.to_numeric(maf['majorAlleleFrequency'], errors='coerce')
    maf['expectedMAF'] = pd.to_numeric(maf['expectedMAF'], errors='coerce')
    maf = maf.dropna(subset=['majorAlleleFrequency', 'expectedMAF']).copy()
    maf['MAF_excess'] = maf['majorAlleleFrequency'] - maf['expectedMAF']
    maf[['tumour_original_source', 'gene', 'cnv_state', 'dm_status', 'aseResults',
         'majorAlleleFrequency', 'expectedMAF', 'MAF_excess']].to_csv(
             str(outdir / 'observed_minus_expected_MAF_gene_sample.tsv'),
             sep='\t', index=False)
    maf.groupby(['cnv_state', 'dm_status', 'aseResults'])['MAF_excess'].agg(
        ['count', 'median', 'mean']).reset_index().to_csv(
            str(outdir / 'observed_minus_expected_MAF_summary.tsv'),
            sep='\t', index=False)
    print('\nCMH analysis summary:\n{}'.format(summary.to_string(index=False)))
    print('\nBalanced-CN sensitivity: {} informative samples with >= {} aDM genes'.format(
        len(sensitivity), args.min_adm))
    print('\nOutput: {}'.format(outdir))
    print('Note: CMH stratifies by sample but does not model dependence from recurrent genes.')


if __name__ == '__main__':
    main()
