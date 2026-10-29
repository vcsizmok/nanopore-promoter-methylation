# Promoter Differential Methylation Analysis

## Overview

This directory contains scripts for identifying and summarizing promoter differential methylation (DM) using haplotype-resolved nanopore methylation data from POG tumour and matched blood samples.

The analysis consists of three stages:

1. Compare CpG-level methylation differences between haplotypes (HP1 and HP2) against cohort-specific reference statistics.
2. Apply DMR-level filtering criteria and annotate tumour DM events using Ploidetect copy-number results.
3. Summarize gene-level DM recurrence and apply additional filters based on normal methylation, lymphoid tissue expression, and HGNC gene classification.

POG and BLOOD samples are analyzed separately, with cohort-specific filtering procedures where appropriate.

## Scripts

| Script | Description |
|---|---|
| `compare_methylation_to_cohort.py` | Compare haplotype methylation differences against cohort-specific CpG statistics |
| `filter_promoter_DM.py` | Apply per-sample promoter DM filters, annotate POG events with Ploidetect, and calculate per-sample DM counts |
| `summarize_promoter_DM.py` | Summarize gene recurrence, exclude hypermethylated POG samples, and apply downstream gene-level filters |

Scripts should be executed in the order listed above.

## Input Data

Input datasets are publicly available at:

https://www.bcgsc.ca/downloads/promoter_methylation_nanopore/

The following datasets are used:

| Dataset | Relative directory or file |
|---|---|
| HP1 methylation frequencies | `methylation_frequency/{POG\|BLOOD}/HP1/` |
| HP2 methylation frequencies | `methylation_frequency/{POG\|BLOOD}/HP2/` |
| Cohort CpG statistics | `methylation_frequency/{POG\|BLOOD}/haplotype_difference_statistics/` |
| DSS promoter DMRs | `dss/{POG\|BLOOD}/filtered_dss_refTSS_hg38_single_best_transcript_intersect_50percent_overlap.tsv.gz` |
| Per-sample cohort comparison | `dmr_comparison_to_cohort/{POG\|BLOOD}/` |
| RefTSS promoter annotation | `annotations/refTSS/refTSS_hg38_single_best_transcript.bed.gz` |
| Imprinted genes | `annotations/imprinted_genes/hg38_mapped_imprinted_genes_geneimprint_db.tsv.gz` |
| Normal methylation classifications | `annotations/normal_methylation/{POG\|BLOOD}/classification_in_normal_for_genes_w_DMR.tsv.gz` |
| Lymphoid-elevated genes | `annotations/lymphoid_genes/tissue_category_rna_lymphoid_elevated_HPA.tsv.gz` |
| HGNC annotations | `annotations/hgnc/hgnc_complete_set.production.tsv.gz` |
| Ploidetect gene-level copy number | `copy_number/POG/<sample_id>_cna_genes.bed.gz` |

All reference files are gzip-compressed.

The Ploidetect and lymphoid gene references are used only for POG tumour samples.

## 1. Comparison to Cohort Methylation Statistics

### Method

For each CpG covered by both haplotypes, the absolute methylation frequency difference is calculated:

\[
\Delta_{\mathrm{CpG}} = |\mathrm{HP1}_{\mathrm{CpG}}-\mathrm{HP2}_{\mathrm{CpG}}|
\]

A CpG is classified as an outlier when:

\[
\Delta_{\mathrm{CpG}} > 2 \times SD_{\mathrm{cohort,CpG}}
\]

The cohort reference statistics are calculated separately for POG and BLOOD samples.

The analysis is restricted to DSS DMRs overlapping at least 50% of an annotated gene promoter.

For each DMR, the proportion of outlier CpGs is calculated:

\[
\mathrm{Proportion\ of\ outlier\ CpGs} =
\frac{\mathrm{Number\ of\ outlier\ CpGs}}
{\mathrm{Number\ of\ evaluable\ CpGs}}
\]

Only CpGs with methylation-frequency measurements for both haplotypes and available cohort reference statistics contribute to the comparison.

### Usage

POG tumour samples:

```bash
python scripts/02_promoter_DM/compare_methylation_to_cohort.py \
  --sample-type POG \
  --data-dir ./data/promoter_methylation_nanopore \
  --output-dir ./results/POG
```

BLOOD samples:

```bash
python scripts/02_promoter_DM/compare_methylation_to_cohort.py \
  --sample-type BLOOD \
  --data-dir ./data/promoter_methylation_nanopore \
  --output-dir ./results/BLOOD
```

### Output

One gzip-compressed TSV file per sample:

```text
<sample_id>_dmr_comparison_to_cohort.tsv.gz
```

| Column | Description |
|---|---|
| `dmr_hugo` | Gene symbol associated with the DMR |
| `meanMethy1_dmr_file` | Mean methylation level for group 1 from DSS |
| `meanMethy2_dmr_file` | Mean methylation level for group 2 from DSS |
| `nCpG_dmr_file` | Number of CpGs reported by DSS |
| `nCpG_methyl_freq_file` | Number of evaluable CpGs |
| `nCpG_outlier` | Number of CpGs with absolute HP1–HP2 difference exceeding 2 × cohort SD |
| `perc_CpG_outlier` | Proportion of evaluable CpGs classified as outliers |

DMR coordinates and the associated gene symbol are retained in the row identifier.

Intermediate chromosome-sorted methylation-frequency files are not required.

## 2. Per-Sample Promoter DM Filtering

Script: `filter_promoter_DM.py`

The following criteria are applied to the DMR-level comparison results.

### Filtering criteria

1. **Autosomal DMRs:** Exclude DMRs on chromosomes X and Y.
2. **CpG outlier proportion:** Require at least 75% of evaluable CpGs within the DMR to be outliers.
3. **CpG representation:** Require at least 75% of CpGs reported in the DSS DMR to be represented in the methylation-frequency data.
4. **Promoter CpG coverage:** Require the number of evaluable CpGs in the DMR to be at least 30% of the CpGs observed within the corresponding promoter window.
5. **Imprinted genes:** Exclude genes present in the imprinted-gene reference.

Promoter windows are defined relative to the annotated transcription start site (TSS), extending 1,500 bp upstream and 1,000 bp downstream, accounting for gene strand.

For promoter CpG coverage, CpGs are counted separately in HP1 and HP2 methylation-frequency files. The larger count is used as the denominator for the 30% requirement.

### Ploidetect annotation (POG only)

For POG tumour samples, filtered promoter DM events are annotated using gene-level Ploidetect copy-number results.

The following fields are retained:

- `zygosity`
- `LOH`
- `A`
- `B`
- `segment_depth`

The `zygosity` field provides the HET/HOM classification used in subsequent gene-level recurrence analysis.

When multiple Ploidetect records are available for the same gene, the first record is retained, consistent with the original analysis.

Ploidetect annotation is not performed for BLOOD samples.

### Per-sample DM counts

For POG, the script also calculates the number of DMR-promoter observations passing the autosomal, CpG outlier, CpG representation, and promoter coverage requirements.

These counts are calculated before imprinted-gene exclusion and before downstream normal methylation, lymphoid, or protein-coding filtering.

They are used to identify hypermethylated tumour samples during gene-level summarization.

### Usage

```bash
python scripts/02_promoter_DM/filter_promoter_DM.py \
  --sample-type POG \
  --data-dir ./data/promoter_methylation_nanopore \
  --output-dir ./results/POG
```

For BLOOD, replace `POG` with `BLOOD` in the sample type and output directory.

## 3. Gene-Level Summarization and Filtering

Script: `summarize_promoter_DM.py`

This script summarizes promoter DM recurrence across samples and applies additional gene-level filtering criteria.

### 3.1 Exclusion of hypermethylated POG samples

Hypermethylated samples are identified from the initial per-sample promoter DM counts.

The threshold is defined as:

\[
T = \operatorname{mean}(N_{\mathrm{DM}}) + 2 \times SD(N_{\mathrm{DM}})
\]

where the standard deviation is calculated with `ddof=0`.

POG samples with counts greater than or equal to this threshold are excluded from gene-level recurrence analysis.

The threshold is calculated once, before normal methylation, lymphoid, and protein-coding gene filtering.

This exclusion is not applied to BLOOD samples.

### 3.2 Gene-level recurrence

DM events are summarized by gene across retained samples.

Genes with symbols containing `-AS1`, `-AS2`, or `-AS3` are excluded.

For POG, recurrence is summarized separately by Ploidetect zygosity:

| Column | Description |
|---|---|
| `gene` | Gene symbol |
| `total_count` | Total number of HET and HOM DM observations |
| `pogs_w_DMs` | Samples contributing DM observations |
| `aDM_count` | Number of HET observations |
| `pogs_w_aDMs` | Samples contributing HET observations |
| `hAM_count` | Number of HOM observations |
| `pogs_w_hAMs` | Samples contributing HOM observations |

The historical `hAM` column names are retained for compatibility with the original analysis.

For BLOOD, the summary contains `gene`, `total_count`, and `pogs_w_DMs`, without separate HET/HOM classification.

### 3.3 Normal methylation filtering

Applied to both POG and BLOOD.

Gene-level DM results are compared against the corresponding normal methylation classification table.

A gene is retained only when:

- It is present in the normal methylation classification table.
- Its `category` is not `methylated_on_one`.

The corresponding normal methylation category is retained in the `in_normal` column.

Genes absent from the classification table are not retained.

### 3.4 Lymphoid gene filtering (POG only)

Genes with elevated expression in lymphoid tissues, as identified by the Human Protein Atlas, are excluded from the POG gene-level summary.

The reference contains genes classified as having elevated RNA expression in lymphoid tissues relative to other tissue types.

This step is not applied to BLOOD samples.

### 3.5 HGNC annotation and protein-coding gene filtering

Applied to both POG and BLOOD.

Remaining genes are annotated using the HGNC production reference.

Genes not found under their original symbols are resolved using the curated symbol mappings from the original analysis. A small number of genes absent from the production reference are manually classified as pseudogenes.

The following columns are added:

| Column | Description |
|---|---|
| `production_hgnc_gene` | HGNC production symbol or manually classified original symbol |
| `locus_group` | HGNC locus classification |

The original gene symbol is preserved in the `gene` column.

Only genes with:

```text
locus_group == "protein-coding gene"
```

are retained in the final protein-coding gene summary.

### Recurrence summaries

At the final filtering stage, the number of genes with `total_count > 1` and `total_count > 10` is also calculated.

These recurrence thresholds are used for summary statistics and do not alter the full protein-coding output.

### Usage

POG:

```bash
python scripts/02_promoter_DM/summarize_promoter_DM.py \
  --sample-type POG \
  --data-dir ./data/promoter_methylation_nanopore \
  --output-dir ./results/POG
```

BLOOD:

```bash
python scripts/02_promoter_DM/summarize_promoter_DM.py \
  --sample-type BLOOD \
  --data-dir ./data/promoter_methylation_nanopore \
  --output-dir ./results/BLOOD
```

## Cohort-Specific Analysis Differences

| Analysis step | POG | BLOOD |
|---|---|---|
| Comparison against cohort CpG statistics | Yes | Yes |
| Per-sample promoter DM filtering | Yes | Yes |
| Ploidetect copy-number annotation | Yes | No |
| Hypermethylated sample exclusion | Yes | No |
| HET/HOM recurrence classification | Yes | No |
| Normal methylation filtering | Yes | Yes |
| Lymphoid-elevated gene exclusion | Yes | No |
| HGNC protein-coding filtering | Yes | Yes |

## Reproducibility Notes

- POG and BLOOD samples use separate cohort methylation reference statistics.
- The CpG outlier threshold is strictly greater than two times the cohort standard deviation, without adding the cohort mean.
- The hypermethylation threshold is calculated before downstream gene-level filtering.
- Original gene symbols are preserved during HGNC harmonization.
- Intermediate results should be retained to allow comparison of gene and sample counts at successive filtering stages.
- Public reference files should be used with the corresponding versions of the analysis scripts.

## Dependencies

- Python 3
- pandas
- numpy