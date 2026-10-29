# Allele-Specific Promoter Methylation in Metastatic Cancers

This repository contains analysis scripts developed for a study investigating allele-specific promoter DNA methylation in metastatic cancers using nanopore long-read sequencing.

The analyses characterize the prevalence and recurrence of allele-specific promoter methylation, its relationship with gene expression and allele-specific expression, and the genomic and epigenomic features of recurrent promoter methylation hotspots.

## Repository structure

| Directory | Description |
|---|---|
| `01_methylation_frequency/` | Analysis of promoter methylation frequencies across tumour samples and cancer types. |
| `02_promoter_DM/` | Identification, characterization, and recurrence analysis of differentially methylated promoters. |
| `03_ASE_analysis/` | Integration of promoter differential methylation with allele-specific expression and copy-number information. |
| `04_hotspot_analysis/` | Identification of genomic regions enriched for promoter differential methylation and characterization of their sequence, chromatin, and regulatory features. |

Each directory contains analysis scripts and a README describing its inputs, methods, dependencies, and execution instructions.

## Analysis overview

**Promoter differential methylation**

Allele-specific differential methylation (aDM) and hemizygous differential methylation (hDM) were characterized across metastatic tumour samples. Analyses include promoter-level recurrence, tumour-type distributions, and comparisons with reference normal methylation datasets.

**Allele-specific expression**

Promoter methylation results were integrated with allele-specific RNA expression and copy-number information to evaluate associations between differential methylation and allelic transcriptional imbalance.

**Promoter methylation hotspots**

A permutation-based approach was used to identify genomic regions enriched for promoters exhibiting differential methylation. The resulting hotspots were characterized using CpG sequence features, ChromHMM chromatin states, CTCF binding, histone modifications, and topologically associating domain (TAD) annotations.

## Software and dependencies

The analyses were developed primarily in Python and use additional bioinformatics tools and reference annotations as described in the README files within each analysis directory.

Dependencies vary by analysis. Please consult the corresponding directory README for specific requirements and execution commands.

## Input data and reference annotations

The analyses use promoter methylation summaries, genomic annotations, and publicly available reference datasets.

Shared genomic annotations are available from:

https://www.bcgsc.ca/downloads/promoter_methylation_nanopore/

Additional reference datasets were obtained from resources including ENCODE and the Roadmap Epigenomics Project.

Large reference datasets, sequencing files, and intermediate analysis outputs are not included in this repository. Scripts provide configurable paths for these external resources.

## Reproducibility

The analysis scripts are organized by analytical objective. Where applicable, permutation-based analyses use documented parameters and random seeds.

For the hotspot functional-characterization analyses, the original set of 44 promoter methylation hotspots is provided as the input for downstream analyses.

See the individual analysis directories for detailed instructions on reproducing results.

## Citation

Publication details will be added when available.