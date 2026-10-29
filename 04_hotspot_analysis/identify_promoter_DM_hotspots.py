#!/usr/bin/env python3
"""Chromosome-stratified promoter-DM hotspot scan (Python 3.7+).

IMPORTANT: Candidate regions are selected and tested using the same observed
labels. Region-level empirical p-values and BH q-values are exploratory and
NOT selection-adjusted. This reproduces the original analysis design.
"""
import argparse
import os
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--promoters', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'input_promoters_coordinates_methylation.tsv'),
                   help='Prepared promoter table (default: 04_hotspot_analysis/data/input_promoters_coordinates_methylation.tsv)')
    p.add_argument('--outdir', default=os.path.join(os.path.dirname(__file__), 'results', 'hotspot_detection'))
    p.add_argument('--window', type=int, default=50000)
    p.add_argument('--step', type=int, default=5000)
    p.add_argument('--min-promoters', type=int, default=3)
    p.add_argument('--window-permutations', type=int, default=5000)
    p.add_argument('--region-permutations', type=int, default=10000)
    p.add_argument('--candidate-regions', default=None,
                   help='Existing candidate-region TSV/BED (e.g. original methylation_regions_from_horspot_analysis.bed). Skips window screening and merging; recomputes p/q values for these exact regions.')
    p.add_argument('--scan-p', type=float, default=0.05)
    p.add_argument('--merge-gap', type=int, default=50000)
    p.add_argument('--seed', type=int, default=12345)
    return p.parse_args()


def make_windows(df, window, step, minimum):
    rows = []
    for chrom, group in df.groupby('chr'):
        positions = np.sort(group['midpoint'].to_numpy(dtype=np.int64))
        for start in range(int(positions.min()), int(positions.max()), step):
            end = start + window
            n = int(np.searchsorted(positions, end, 'left') - np.searchsorted(positions, start, 'left'))
            if n >= minimum:
                rows.append((chrom, start, end, n))
    return pd.DataFrame(rows, columns=['chr', 'start', 'end', 'n_promoters'])


def count_observed(intervals, promoters, inclusive=False):
    counts = np.zeros(len(intervals), dtype=np.int64)
    for chrom, group in intervals.groupby('chr', sort=False):
        subset = promoters.loc[promoters['chr'] == chrom].sort_values('midpoint')
        positions = subset['midpoint'].to_numpy(dtype=np.int64)
        values = subset['methylated'].to_numpy(dtype=np.int64)
        left = np.searchsorted(positions, group['start'].to_numpy(dtype=np.int64), side='left')
        right = np.searchsorted(positions, group['end'].to_numpy(dtype=np.int64), side='right' if inclusive else 'left')
        csum = np.concatenate(([0], np.cumsum(values)))
        counts[group.index.to_numpy()] = csum[right] - csum[left]
    return counts


def perm_pvalues(intervals, promoters, n_perm, rng, inclusive=False):
    observed = count_observed(intervals, promoters, inclusive)
    exceed = np.zeros(len(intervals), dtype=np.int64)
    chrom_data = {}
    for chrom, group in promoters.groupby('chr'):
        sorted_group = group.sort_values('midpoint')
        chrom_data[chrom] = (sorted_group['midpoint'].to_numpy(dtype=np.int64), sorted_group['methylated'].to_numpy(dtype=np.int64))
    lookups = []
    for chrom, group in intervals.groupby('chr', sort=False):
        positions, values = chrom_data[chrom]
        left = np.searchsorted(positions, group['start'].to_numpy(dtype=np.int64), 'left')
        right = np.searchsorted(positions, group['end'].to_numpy(dtype=np.int64), 'right' if inclusive else 'left')
        lookups.append((group.index.to_numpy(), chrom, left, right))
    for i in range(n_perm):
        shuffled = {chrom: np.concatenate(([0], np.cumsum(rng.permutation(values)))) for chrom, (positions, values) in chrom_data.items()}
        for indices, chrom, left, right in lookups:
            csum = shuffled[chrom]
            exceed[indices] += (csum[right] - csum[left] >= observed[indices])
        if (i + 1) % 1000 == 0:
            print('Permutations: {}/{}'.format(i + 1, n_perm))
    return observed, (exceed + 1.0) / (n_perm + 1.0)


def merge_regions(windows, gap):
    regions = []
    for chrom, group in windows.groupby('chr'):
        current = None
        for row in group.sort_values('start').itertuples(index=False):
            if current is None:
                current = [chrom, int(row.start), int(row.end)]
            elif int(row.start) <= current[2] + gap:
                current[2] = max(current[2], int(row.end))
            else:
                regions.append(current)
                current = [chrom, int(row.start), int(row.end)]
        if current is not None:
            regions.append(current)
    return pd.DataFrame(regions, columns=['chr', 'start', 'end'])


def main():
    a = parse_args()
    if a.window <= 0 or a.step <= 0 or a.window_permutations <= 0 or a.region_permutations <= 0:
        raise ValueError('Window, step and permutation counts must be positive')
    os.makedirs(a.outdir, exist_ok=True)
    if not os.path.isfile(a.promoters):
        raise IOError('Promoter input not found: {}. Run prepare_hotspot_promoter_input.py first.'.format(a.promoters))
    print('Promoter input: {}'.format(a.promoters))
    df = pd.read_csv(a.promoters, sep='\t')
    if 'symbol' not in df.columns and df.index.name == 'symbol':
        df = df.reset_index()
    required = ['chrom', 'promoter_start', 'promoter_end', 'methylated']
    missing = [x for x in required if x not in df.columns]
    if missing:
        raise ValueError('Missing columns: ' + ', '.join(missing))
    df = df.rename(columns={'chrom': 'chr'}).copy()
    df['midpoint'] = ((pd.to_numeric(df['promoter_start']) + pd.to_numeric(df['promoter_end'])) // 2).astype(np.int64)
    df['methylated'] = pd.to_numeric(df['methylated'], errors='raise')
    if df['methylated'].isna().any() or not df['methylated'].isin([0, 1]).all():
        raise ValueError('methylated must contain only binary 0/1 labels')
    rng = np.random.RandomState(a.seed)
    if a.candidate_regions is not None:
        if not os.path.isfile(a.candidate_regions):
            raise IOError('Candidate-region file not found: {}'.format(a.candidate_regions))
        regions = pd.read_csv(a.candidate_regions, sep='\t')
        if not set(['chr', 'start', 'end']).issubset(regions.columns):
            raise ValueError('Candidate regions must contain chr, start, end columns')
        regions = regions[['chr', 'start', 'end']].copy()
        if regions.empty or regions[['chr', 'start', 'end']].isna().any().any():
            raise ValueError('Candidate regions are empty or contain missing coordinates')
        regions[['start', 'end']] = regions[['start', 'end']].astype(np.int64)
        if (regions['end'] < regions['start']).any():
            raise ValueError('Candidate region end must be >= start')
        if regions.duplicated(['chr', 'start', 'end']).any():
            raise ValueError('Duplicate candidate regions')
        print('Reusing {} candidate regions from {}'.format(len(regions), a.candidate_regions))
    else:
        windows = make_windows(df, a.window, a.step, a.min_promoters)
        if windows.empty:
            raise ValueError('No windows meet the minimum promoter count')
        print('Testing {} eligible windows'.format(len(windows)))
        windows['meth_count'], windows['pval'] = perm_pvalues(windows, df, a.window_permutations, rng, inclusive=False)
        windows.to_csv(os.path.join(a.outdir, 'window_screening_results.tsv'), sep='\t', index=False)
        regions = merge_regions(windows.loc[windows['pval'] < a.scan_p], a.merge_gap)
        if regions.empty:
            print('No candidate regions passed nominal window screening')
            pd.DataFrame(columns=['chr', 'start', 'end', 'meth_count', 'n_promoters', 'pval', 'qval']).to_csv(os.path.join(a.outdir, 'candidate_regions.tsv'), sep='\t', index=False)
            return
    # Original analysis uses inclusive ends for region counts, unlike screening windows.
    for chrom, group in regions.groupby('chr'):
        positions = np.sort(df.loc[df['chr'] == chrom, 'midpoint'].to_numpy(dtype=np.int64))
        regions.loc[group.index, 'n_promoters'] = np.searchsorted(positions, group['end'].to_numpy(dtype=np.int64), 'right') - np.searchsorted(positions, group['start'].to_numpy(dtype=np.int64), 'left')
    regions = regions.loc[regions['n_promoters'] >= a.min_promoters].reset_index(drop=True)
    if regions.empty:
        print('No merged regions meet minimum promoter count')
        return
    print('Testing {} candidate regions with {} permutations (seed={})'.format(len(regions), a.region_permutations, a.seed))
    regions['meth_count'], regions['pval'] = perm_pvalues(regions, df, a.region_permutations, rng, inclusive=True)
    regions['qval'] = multipletests(regions['pval'], method='fdr_bh')[1]
    regions.to_csv(os.path.join(a.outdir, 'candidate_regions.tsv'), sep='\t', index=False)
    hotspots = regions.loc[regions['qval'] < 0.05].copy()
    hotspots['meth_fraction'] = hotspots['meth_count'] / hotspots['n_promoters']
    hotspots['hotspot_score'] = -np.log10(hotspots['qval'])
    hotspots[['chr', 'start', 'end', 'meth_count', 'n_promoters', 'meth_fraction', 'pval', 'qval', 'hotspot_score']].to_csv(os.path.join(a.outdir, 'methylated_promoter_hotspots.bed'), sep='\t', index=False)
    print('{} exploratory hotspots; outputs: {}'.format(len(hotspots), a.outdir))


if __name__ == '__main__':
    main()
