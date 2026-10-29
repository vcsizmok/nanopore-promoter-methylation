# Allele-specific expression (ASE) and promoter differential methylation

This directory contains the scripts used to integrate allele-specific expression (ASE) results from IMPALA with promoter differential methylation (DM) and copy-number annotations, and to assess associations between ASE, promoter DM, and copy-number state.

The workflow uses **IMPALA's original ASE/BAE classifications**, without imposing an additional adjusted-*p*-value threshold. It is designed for **Python 3.7** and uses the validated protein-coding promoter DM results from `02_promoter_DM`.

## Directory layout

```text
03_ASE_analysis/
├── data/
│   ├── summary_protein_coding.tsv
│   └── hypermethylated_samples.tsv
├── integrate_IMPALA_ASE_DM_CNV.py
├── compare_ASE_BAE_frequencies.py
├── analyze_directional_aDM_ASE.py
├── analyze_aDM_ASE_by_copy_number.py
├── README.md
└── results/                         # generated when scripts are run
```

The two files in `data/` are derived from the validated `02_promoter_DM` analysis. They are included as explicit, fixed inputs to this workflow. `summary_protein_coding.tsv` provides gene-level lists of samples with aDM and hDM; `hypermethylated_samples.tsv` contains the `pogid` column identifying the 15 samples excluded from promoter DM recurrence analysis.

## Public IMPALA input data

The cohort-selected IMPALA summary tables are deposited at:

https://www.bcgsc.ca/downloads/promoter_methylation_nanopore/ase/IMPALA/

Expected directory layout:

```text
ase/IMPALA/
├── <sample_1>/summaryTable.tsv.gz
├── <sample_2>/summaryTable.tsv.gz
└── ...
```

Each table must contain the IMPALA fields required by the integration script, including `gene`, `aseResults`, `allele1IsMajor`, `majorAlleleFrequency`, `cnv_state`, `expectedMAF`, and `methyl_state`. The integration script also accepts uncompressed `summaryTable.tsv` files.

**Cohort definition:** The deposited IMPALA directories have already been restricted to the study cohort. A separate cohort sample table is therefore not required by this workflow. The integration script then excludes the same **15 hypermethylated samples** used by `02_promoter_DM`, including `POG1098-OCT-1`. The number of samples with usable IMPALA data may be smaller than the 363 samples retained in the promoter DM recurrence analysis.

## Software requirements

- Python 3.7
- NumPy
- pandas
- SciPy
- matplotlib
- statsmodels

Install dependencies in an environment appropriate for your system. For reproducibility, record the actual package versions used for the analysis; no specific versions are asserted here.

## Run the analysis

Run the commands below from `03_ASE_analysis/`. Default input/output paths are resolved relative to each script's location, except for the IMPALA input directory, whose default points to the GSC downloads filesystem.

### 1. Integrate IMPALA, promoter DM, and copy number

```bash
python integrate_IMPALA_ASE_DM_CNV.py
```

To use a downloaded copy of the IMPALA data instead of the GSC filesystem:

```bash
python integrate_IMPALA_ASE_DM_CNV.py \
  --ase-dir /path/to/ase/IMPALA \
  --dm-file data/summary_protein_coding.tsv \
  --excluded-samples data/hypermethylated_samples.tsv \
  --output-dir results
```

The integration script:

- Retains only records classified by IMPALA as `ASE` or `BAE`, without additional `padj` filtering.
- Uses `pogs_w_aDMs` and `pogs_w_hAMs` from the validated DM table to assign sample-specific `aDM`, `hDM`, `both_aDM_hDM`, or `no_DM` status.
- Removes all records for the 15 hypermethylated samples specified by `pogid`.
- Preserves `balance`, `imbalance`, and `LOH` copy-number states; unrecognized or unavailable states remain missing, **not balanced**.
- Audits duplicate gene records within each IMPALA sample. All rows for duplicated genes are excluded from the primary integrated table and retained in a separate audit file.

Main outputs under `results/`:

| File | Description |
| --- | --- |
| `IMPALA_ASE_DM_CNV_gene_sample_table.tsv` | One retained observation per gene and sample, with IMPALA measurements, DM status, and CN state |
| `IMPALA_ASE_BAE_sample_summary.tsv` | ASE/BAE counts by DM and CN category for each sample |
| `ASE_BAE_cnvState_count/<sample>_ASE_BAE_cnvState_count.tsv` | Individual sample count tables |
| `IMPALA_duplicate_gene_records.tsv` | Duplicate IMPALA gene records excluded from the primary analysis |
| `IMPALA_processing_summary.tsv` | Processing status and observation counts |

### 2. Compare ASE and BAE frequencies

```bash
python compare_ASE_BAE_frequencies.py
```

The script compares, within each sample, the percentages of ASE versus BAE genes that exhibit **aDM, LOH, imbalanced CN, or balanced CN**. The primary statistical comparisons are **paired Wilcoxon signed-rank tests** (`zero_method='zsplit'`). Samples without usable observations for either expression class are excluded pairwise.

For aDM, the denominator is genes with recognized DM status; for CN analyses, the denominator is genes with recognized CN state. These denominators are calculated separately. Genes classified as `both_aDM_hDM` are included in the aDM category for this comparison.

A pooled Fisher exact test for aDM is also reported **descriptively**; it does not account for clustering within samples.

Outputs are written to `results/ASE_BAE_cnvState_count/wilcoxon_unfiltered/`, including `ASE_BAE_paired_Wilcoxon_results.tsv`, per-sample counts and percentages, and PNG/PDF plots.

### 3. Test directional concordance between aDM and ASE

```bash
python analyze_directional_aDM_ASE.py
```

The analysis selects observations classified as both `aDM` and `ASE`. IMPALA's `methyl_state` is interpreted as **HP1 minus HP2 methylation**, and `allele1IsMajor` indicates whether **HP1 is more highly expressed**. An observation is directionally concordant when the less-methylated haplotype is more highly expressed.

The **primary analysis** uses genes in **balanced-CN regions**. Secondary analyses include all informative CN states, imbalanced CN, and LOH. Sensitivity analyses restrict balanced-CN observations to an absolute methylation difference of at least **0.20** or samples with at least **five** informative observations. Both thresholds can be changed via `--min-abs-methyl` and `--min-sample-events`.

The script reports exact binomial tests with Clopper–Pearson 95% confidence intervals and sample-level sign tests. **Gene–sample-level binomial tests assume independent observations** and do not account for within-sample or recurrent-gene dependence; interpret them alongside the sample-level summaries.

Outputs are written to `results/directional_aDM_ASE_concordance/`, including `directional_concordance_summary.tsv`, `per_sample_directional_concordance.tsv`, `sample_level_sign_tests.tsv`, data-availability audits, and PNG/PDF figures.

### 4. Test ASE enrichment among aDM genes by copy-number state

```bash
python analyze_aDM_ASE_by_copy_number.py
```

The primary comparison is **aDM versus no promoter DM in balanced-CN regions**, using a **Cochran–Mantel–Haenszel (CMH) test stratified by sample**. Secondary analyses use imbalanced CN, LOH, and all known CN states. Only observations classified as `aDM` or `no_DM` are included; hDM and combined aDM/hDM observations are not included in this analysis.

A sample contributes to a CMH test only when both DM categories and both ASE/BAE categories are represented. Sample-specific odds ratios use a **0.5 Haldane–Anscombe correction**, but the CMH tables use the uncorrected counts. The analysis also reports sample-level sign tests and, where estimable, heterogeneity tests.

The sensitivity analysis retains samples with at least **five aDM genes within balanced-CN regions** (`--min-adm 5` by default). The script additionally summarizes observed minus expected major allele frequency (`majorAlleleFrequency - expectedMAF`).

Outputs are written to `results/copy_number_stratified_aDM_ASE/`, including:

- `CMH_analysis_summary.tsv`
- `CMH_analysis_summary_with_sensitivity.tsv`
- `balanced_CNV_aDM_ge5_CMH_summary.tsv`
- `sample_specific_aDM_ASE_results.tsv`
- `balanced_CNV_aDM_ge5_sample_results.tsv`
- `observed_minus_expected_MAF_gene_sample.tsv`
- `observed_minus_expected_MAF_summary.tsv`
- PNG/PDF figures for sample-level odds ratios, paired ASE frequencies, and CMH effect sizes

The CMH analysis accounts for stratification by tumour sample but does not explicitly model dependence due to recurrent genes.

## Interpretation and reproducibility notes

1. **The promoter DM calls are not recalculated here.** Sample-specific aDM/hDM assignments come from `02_promoter_DM`'s validated protein-coding summary.
2. **Manual IGV-reviewed BRCA1/RAD51C observations** that did not meet the automated promoter DM thresholds are not added to the automated ASE integration. These cases should be described separately in the manuscript, where applicable.
3. **The same 15-sample exclusion set** is applied at integration. Downstream statistical scripts use the already-filtered integrated table and do not exclude additional samples.
4. **Missing CN is not balanced CN.** Analyses restricted to known CN states exclude unavailable CN annotations from their denominators.
5. **IMPALA classifications are retained as supplied.** No additional `padj` cutoff is applied by these scripts.
6. **Counts represent retained gene–sample observations**, following the duplicate-record exclusion policy. These are not necessarily counts of unique genes across the cohort.
7. For publication-grade reproducibility, retain the exact versions of the deposited IMPALA tables, the two DM input TSVs, and the Python dependencies used to generate the reported results.

Each script supports `--help` to display its configurable inputs and output paths.
