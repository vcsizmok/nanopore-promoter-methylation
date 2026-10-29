#!/usr/bin/env python3

"""
Calculate the mean and standard deviation of absolute haplotype
methylation-frequency differences at each CpG across samples.

Inputs:
    HP1 and HP2 methylation frequency files intersected with
    promoter-overlapping DSS regions.

Outputs:
    One TSV per chromosome containing:
        CpG site
        mean_methyl_freq
        std_methyl_freq

Standard deviation uses ddof=0, matching the original code.
"""

import argparse
import gzip
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd


CHROMOSOMES = [
    f"chr{i}" for i in range(1, 23)
] + ["chrX", "chrY"]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate cohort-level HP1/HP2 methylation differences"
    )
    parser.add_argument(
        "--sample-type",
        choices=["POG", "BLOOD"],
        required=True,
        help="Sample cohort to process"
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing POG/ and BLOOD/ subdirectories"
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory for output files"
    )
    return parser.parse_args()


def read_methylation_file(filename, chromosome):
    """Read CpG methylation frequencies for one chromosome."""
    cpg_freq = {}

    opener = gzip.open if str(filename).endswith(".gz") else open

    with opener(filename, "rt") as handle:
        for line in handle:
            columns = line.rstrip("\n").split("\t")

            if columns[0] != chromosome:
                continue

            cpg = f"{columns[0]}_{columns[1]}_{columns[2]}"
            cpg_freq[cpg] = float(columns[5])

    return cpg_freq


def get_sample_files(input_dir, sample_type):
    cohort_dir = Path(input_dir) / sample_type

    hp1_dir = cohort_dir / "HP1"
    hp2_dir = cohort_dir / "HP2"

    suffix_hp1 = "_HP1_methyl_freq_intersect_w_filtered_dss_regions.tsv.gz"
    suffix_hp2 = "_HP2_methyl_freq_intersect_w_filtered_dss_regions.tsv.gz"

    hp1_files = {
        p.name[:-len(suffix_hp1)]: p
        for p in hp1_dir.glob(f"*{suffix_hp1}")
    }

    hp2_files = {
        p.name[:-len(suffix_hp2)]: p
        for p in hp2_dir.glob(f"*{suffix_hp2}")
    }

    samples = sorted(set(hp1_files) & set(hp2_files))

    print(f"HP1 files: {len(hp1_files)}")
    print(f"HP2 files: {len(hp2_files)}")
    print(f"Paired samples: {len(samples)}")

    missing_hp1 = sorted(set(hp2_files) - set(hp1_files))
    missing_hp2 = sorted(set(hp1_files) - set(hp2_files))

    if missing_hp1:
        print(f"WARNING: {len(missing_hp1)} samples missing HP1")
    if missing_hp2:
        print(f"WARNING: {len(missing_hp2)} samples missing HP2")

    return [
        (sample, hp1_files[sample], hp2_files[sample])
        for sample in samples
    ]


def main():
    args = parse_args()

    output_dir = Path(args.output_dir) / args.sample_type
    output_dir.mkdir(parents=True, exist_ok=True)

    sample_files = get_sample_files(
        args.input_dir,
        args.sample_type
    )

    if not sample_files:
        raise RuntimeError("No paired HP1/HP2 files found")

    for chromosome in CHROMOSOMES:
        print(f"\nProcessing {chromosome}")

        hp_diff_dict = defaultdict(list)

        for sample, hp1_file, hp2_file in sample_files:
            hp1 = read_methylation_file(hp1_file, chromosome)
            hp2 = read_methylation_file(hp2_file, chromosome)

            common_cpgs = hp1.keys() & hp2.keys()

            for cpg in common_cpgs:
                diff = abs(hp1[cpg] - hp2[cpg])
                hp_diff_dict[cpg].append(diff)

        results = {
            cpg: [
                np.mean(differences),
                np.std(differences, ddof=0)
            ]
            for cpg, differences in hp_diff_dict.items()
        }

        df = pd.DataFrame.from_dict(
            results,
            orient="index",
            columns=["mean_methyl_freq", "std_methyl_freq"]
        )

        output_file = (
            output_dir /
            f"{chromosome}_mean_std_methyl_freq_diff_cpgs.tsv"
        )

        df.to_csv(output_file, sep="\t", index=True)

        print(f"Saved {len(df)} CpGs to {output_file}")


if __name__ == "__main__":
    main()
