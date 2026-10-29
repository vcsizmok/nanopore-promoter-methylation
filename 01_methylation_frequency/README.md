# Haplotype Methylation Frequency Differences

## Description

This script calculates the mean and standard deviation of absolute methylation frequency differences between haplotypes (HP1 and HP2) at individual CpG sites across samples.

The analysis is performed separately for tumour (POG) and blood (BLOOD) samples. Only CpGs with methylation frequency measurements available for both haplotypes are included.

The standard deviation is calculated using `ddof=0`.

## Input data

Haplotype-specific CpG methylation frequency files intersected with DSS regions covering at least 50% of the promoter window.

Input files are available at:

https://www.bcgsc.ca/downloads/promoter_methylation_nanopore/methylation_frequency/

The expected input directory structure is:

```text
methylation_frequency/
├── POG/
│   ├── HP1/
│   └── HP2/
└── BLOOD/
    ├── HP1/
    └── HP2/
```

## Usage

Run the script from the repository root.

### Tumour samples (POG)

```bash
python scripts/01_methylation_frequency/calculate_haplotype_methylation_difference.py \
  --sample-type POG \
  --input-dir ./data/methylation_frequency \
  --output-dir ./results
```

### Blood samples (BLOOD)

```bash
python scripts/01_methylation_frequency/calculate_haplotype_methylation_difference.py \
  --sample-type BLOOD \
  --input-dir ./data/methylation_frequency \
  --output-dir ./results
```

Download the input files from the BCGSC server and preserve the expected directory structure before running the script.

## Output

One tab-separated file per chromosome containing:

- CpG coordinate (`chr_start_end`)
- `mean_methyl_freq`: Mean absolute HP1–HP2 methylation frequency difference across samples
- `std_methyl_freq`: Standard deviation of the absolute difference across samples

Output files are written separately for POG and BLOOD under the specified output directory.

These chromosome-level statistics are used as reference inputs for downstream promoter differential methylation analyses.