#!/usr/bin/env python3
"""CpG/GC/CpG-island feature analysis using promoter-matched random intervals.

Random regions match chromosome, inclusive interval length, and the exact
number of evaluable promoter midpoints. Python 3.7+; requires bedtools.
"""
import argparse
import os
import shutil
import subprocess
import numpy as np
import pandas as pd

METRICS = ['median_CpG_per_kb', 'median_CpG_observed_expected', 'median_GC_fraction',
           'fraction_CpG_island_overlap', 'median_CpG_island_fraction',
           'median_CpG_island_count', 'median_N_fraction']


def args_parse():
    p = argparse.ArgumentParser(description=__doc__)
    base = os.path.dirname(__file__)
    p.add_argument('--hotspots', default=os.path.join(base, 'data', 'methylated_promoter_hotspots.bed'))
    p.add_argument('--promoters', default=os.path.join(base, 'data', 'input_promoters_coordinates_methylation.tsv'))
    p.add_argument('--chrom-sizes', default='/gsc/www/bcgsc.ca/downloads/promoter_methylation_nanopore/annotations/genome/hg38.chrom.sizes', help='Chromosome sizes TSV, hg38')
    p.add_argument('--fasta', default='/gsc/resources/Homo_sapiens_genomes/hg38_no_alt/genome/fasta/hg38_no_alt.fa', help='Original GSC hg38_no_alt FASTA; override for external runs')
    p.add_argument('--cpg-islands', default='/gsc/www/bcgsc.ca/downloads/promoter_methylation_nanopore/annotations/cpg_islands/cpgIslandExt.hg38.full.bed', help='UCSC hg38 cpgIslandExt BED')
    p.add_argument('--outdir', default=os.path.join(base, 'results', 'CpG_sequence_analysis'))
    p.add_argument('--permutations', type=int, default=10000)
    p.add_argument('--expected-hotspots', type=int, default=44,
                   help='Expected manuscript hotspot count; 0 disables check')
    p.add_argument('--seed', type=int, default=12345)
    return p.parse_args()


def valid_ranges(positions, length, k, chrom_size):
    max_start = chrom_size - length
    if max_start < 0 or k > len(positions):
        return []
    ranges = []
    if k == 0:
        previous = -1
        for position in positions:
            lo, hi = max(0, previous + 1), min(max_start, int(position) - length)
            if lo <= hi:
                ranges.append((lo, hi))
            previous = int(position)
        if previous + 1 <= max_start:
            ranges.append((previous + 1, max_start))
        return ranges
    for i in range(len(positions) - k + 1):
        j = i + k - 1
        first, last = int(positions[i]), int(positions[j])
        if last - first + 1 > length:
            continue
        lo = max(0, last - length + 1)
        hi = min(max_start, first)
        if i > 0:
            lo = max(lo, int(positions[i - 1]) + 1)
        if j + 1 < len(positions):
            hi = min(hi, int(positions[j + 1]) - length)
        if lo <= hi:
            ranges.append((lo, hi))
    return ranges


def sample_start(info, rng):
    index = int(rng.randint(0, info['total']))
    r = int(np.searchsorted(info['cumulative'], index, side='right'))
    previous = 0 if r == 0 else int(info['cumulative'][r - 1])
    return int(info['ranges'][r][0] + index - previous)


def count_promoters(positions, start, end):
    return int(np.searchsorted(positions, end, 'right') - np.searchsorted(positions, start, 'left'))


def write_and_run(command, outfile):
    with open(outfile, 'w') as out:
        subprocess.check_call(command, stdout=out)


def sequence_metrics(filename):
    rows = []
    with open(filename) as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 2:
                raise ValueError('Unexpected getfasta output: ' + line[:100])
            name, seq = fields
            seq = seq.upper()
            na, nc, ng, nt = (seq.count(x) for x in 'ACGT')
            valid = na + nc + ng + nt
            cpg = seq.count('CG')
            rows.append(dict(region_id=name, sequence_length=len(seq), ACGT_length=valid,
                             N_fraction=(len(seq)-valid)/float(len(seq)) if seq else np.nan,
                             n_A=na, n_C=nc, n_G=ng, n_T=nt, n_CpG=cpg,
                             GC_fraction=(nc+ng)/float(valid) if valid else np.nan,
                             CpG_per_kb=1000*cpg/float(valid) if valid else np.nan,
                             CpG_observed_expected=cpg*valid/float(nc*ng) if nc and ng else np.nan))
    return pd.DataFrame(rows)


def set_stats(df):
    return dict(median_CpG_per_kb=np.nanmedian(df['CpG_per_kb']),
                median_CpG_observed_expected=np.nanmedian(df['CpG_observed_expected']),
                median_GC_fraction=np.nanmedian(df['GC_fraction']),
                fraction_CpG_island_overlap=np.mean(df['CpG_island_overlap']),
                median_CpG_island_fraction=np.nanmedian(df['CpG_island_fraction']),
                median_CpG_island_count=np.nanmedian(df['CpG_island_count']),
                median_N_fraction=np.nanmedian(df['N_fraction']))


def summarize(observed, permutations):
    rows = []
    for metric in METRICS:
        obs = float(observed[metric])
        vals = permutations[metric].to_numpy(dtype=float)
        vals = vals[np.isfinite(vals)]
        if not len(vals):
            raise ValueError('No finite permutation values for ' + metric)
        lo = (np.count_nonzero(vals <= obs) + 1.0)/(len(vals)+1)
        hi = (np.count_nonzero(vals >= obs) + 1.0)/(len(vals)+1)
        sd = np.std(vals, ddof=1) if len(vals) > 1 else np.nan
        median = np.median(vals)
        rows.append(dict(metric=metric, observed=obs, random_mean=np.mean(vals), random_median=median,
                         random_2_5_percentile=np.percentile(vals, 2.5), random_97_5_percentile=np.percentile(vals, 97.5),
                         observed_to_random_median_ratio=obs/median if median else np.nan,
                         permutation_z=(obs-np.mean(vals))/sd if sd > 0 else np.nan,
                         empirical_p_enrichment=hi, empirical_p_depletion=lo,
                         empirical_p_two_sided=min(1., 2*min(lo, hi))))
    return pd.DataFrame(rows)


def main():
    a = args_parse()
    if a.permutations < 1:
        raise ValueError('--permutations must be positive')
    for filename in (a.hotspots, a.promoters, a.chrom_sizes, a.fasta, a.cpg_islands):
        if not os.path.isfile(filename):
            raise IOError('Missing input: ' + filename)
    if not shutil.which('bedtools'):
        raise RuntimeError('bedtools not found in PATH')
    os.makedirs(a.outdir, exist_ok=True)
    hotspots = pd.read_csv(a.hotspots, sep='\t').rename(columns={'chr':'chrom'}).reset_index(drop=True)
    promoters = pd.read_csv(a.promoters, sep='\t')
    if 'chrom' not in promoters and 'chr' in promoters:
        promoters = promoters.rename(columns={'chr':'chrom'})
    for col in ('chrom', 'start', 'end', 'n_promoters'):
        if col not in hotspots:
            raise ValueError('Missing hotspot column ' + col)
    for col in ('chrom', 'promoter_start', 'promoter_end'):
        if col not in promoters:
            raise ValueError('Missing promoter column ' + col)
    if hotspots.empty:
        raise ValueError('No hotspots to analyze')
    if a.expected_hotspots and len(hotspots) != a.expected_hotspots:
        raise ValueError('Expected {} reference hotspots, found {} in {}. Use --expected-hotspots 0 only if intentional.'.format(
            a.expected_hotspots, len(hotspots), a.hotspots))
    print('Analyzing {} hotspots from {}'.format(len(hotspots), a.hotspots))
    hotspots[['start','end','n_promoters']] = hotspots[['start','end','n_promoters']].astype(np.int64)
    if (hotspots['end'] < hotspots['start']).any():
        raise ValueError('Invalid hotspot intervals')
    promoters['midpoint'] = ((promoters['promoter_start'].astype(np.int64)+promoters['promoter_end'].astype(np.int64))//2)
    lookup = {chrom:np.sort(group['midpoint'].to_numpy(dtype=np.int64)) for chrom,group in promoters.groupby('chrom')}
    chrom_sizes = pd.read_csv(a.chrom_sizes, sep='\t', header=None, names=['chrom','size'])
    sizes = dict(zip(chrom_sizes['chrom'], chrom_sizes['size'].astype(np.int64)))
    hotspots['inclusive_length'] = hotspots['end']-hotspots['start']+1
    for row in hotspots.itertuples():
        if row.chrom not in lookup or row.chrom not in sizes:
            raise ValueError('Chromosome absent from promoter/chromosome sizes: ' + row.chrom)
        if count_promoters(lookup[row.chrom], row.start, row.end) != row.n_promoters:
            raise ValueError('Observed hotspot promoter-count mismatch at row {}'.format(row.Index))
        if row.end >= sizes[row.chrom]:
            raise ValueError('Hotspot exceeds chromosome boundary')
    cache = {}
    for row in hotspots.itertuples():
        key = (row.chrom, int(row.inclusive_length), int(row.n_promoters))
        if key in cache:
            continue
        ranges = valid_ranges(lookup[row.chrom], key[1], key[2], int(sizes[row.chrom]))
        if not ranges:
            raise ValueError('No valid matched random intervals for {}'.format(key))
        counts = np.asarray([hi-lo+1 for lo,hi in ranges], dtype=np.int64)
        cache[key] = dict(ranges=ranges, cumulative=np.cumsum(counts), total=int(counts.sum()))
        print('{}: {:,} valid start positions'.format(key, cache[key]['total']))
    rng = np.random.RandomState(a.seed)
    rows = []
    for perm in range(a.permutations+1):
        for i, row in enumerate(hotspots.itertuples()):
            length = int(row.inclusive_length)
            start = int(row.start) if perm == 0 else sample_start(cache[(row.chrom,length,int(row.n_promoters))], rng)
            end = start+length-1
            if count_promoters(lookup[row.chrom],start,end) != row.n_promoters:
                raise RuntimeError('Random interval promoter-count mismatch')
            rows.append(dict(chrom=row.chrom,start_inclusive=start,end_inclusive=end,inclusive_length=length,
                             bed_start=start,bed_end=end+1,region_id=('OBS_HS_{:03d}'.format(i+1) if perm == 0 else 'PERM_{:05d}_HS_{:03d}'.format(perm,i+1)),
                             set_type='observed' if perm == 0 else 'random',permutation=perm,hotspot_index=i+1,
                             target_n_promoters=int(row.n_promoters),promoter_count_recalculated=int(row.n_promoters)))
        if perm and perm % 1000 == 0:
            print('Generated {}/{} random sets'.format(perm,a.permutations))
    intervals = pd.DataFrame(rows)
    intervals.to_csv(os.path.join(a.outdir,'all_observed_and_random_regions.tsv'),sep='\t',index=False)
    bedfile = os.path.join(a.outdir,'all_observed_and_random_regions.bed')
    intervals[['chrom','bed_start','bed_end','region_id']].to_csv(bedfile,sep='\t',header=False,index=False)
    seqfile = os.path.join(a.outdir,'all_region_sequences.tsv')
    write_and_run(['bedtools','getfasta','-fi',a.fasta,'-bed',bedfile,'-nameOnly','-tab'],seqfile)
    seq = sequence_metrics(seqfile)
    coveragefile = os.path.join(a.outdir,'all_regions_CpG_island_coverage.tsv')
    write_and_run(['bedtools','coverage','-a',bedfile,'-b',a.cpg_islands],coveragefile)
    coverage = pd.read_csv(coveragefile,sep='\t',header=None)
    if coverage.shape[1] != 8:
        raise ValueError('Expected eight columns from bedtools coverage')
    coverage.columns = ['chrom','bed_start','bed_end','region_id','CpG_island_count','CpG_island_bases','bed_interval_length','CpG_island_fraction']
    coverage['CpG_island_overlap'] = coverage['CpG_island_count'] > 0
    metrics = intervals.merge(seq,on='region_id',how='left',validate='one_to_one').merge(coverage[['region_id','CpG_island_count','CpG_island_bases','bed_interval_length','CpG_island_fraction','CpG_island_overlap']],on='region_id',how='left',validate='one_to_one')
    if len(metrics) != len(intervals) or metrics['sequence_length'].isna().any() or not (metrics['sequence_length'] == metrics['inclusive_length']).all() or not (metrics['bed_interval_length'] == metrics['inclusive_length']).all():
        raise RuntimeError('Sequence or CpG-island coverage length validation failed')
    metrics.to_csv(os.path.join(a.outdir,'all_region_CpG_sequence_metrics.tsv'),sep='\t',index=False)
    observed = metrics.loc[metrics['set_type']=='observed'].copy()
    observed.to_csv(os.path.join(a.outdir,'observed_hotspot_CpG_sequence_metrics.tsv'),sep='\t',index=False)
    observed_stats = set_stats(observed)
    permutation_stats = pd.DataFrame([dict(set_stats(group),permutation=int(perm)) for perm,group in metrics.loc[metrics['set_type']=='random'].groupby('permutation',sort=True)])
    permutation_stats.to_csv(os.path.join(a.outdir,'CpG_sequence_permutation_distributions.tsv'),sep='\t',index=False)
    summary = summarize(observed_stats,permutation_stats)
    summary.to_csv(os.path.join(a.outdir,'CpG_sequence_permutation_summary.tsv'),sep='\t',index=False)
    print(summary.to_string(index=False))
    print('Results: ' + a.outdir)


if __name__ == '__main__':
    main()
