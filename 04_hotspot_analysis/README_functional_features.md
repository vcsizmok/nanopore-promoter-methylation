# Hotspot functional-characterization workflow (corrected paths)

The **original 44 manuscript hotspots** are authoritative; do not regenerate them for this figure. All scripts use paths relative to their own location, so they can be invoked from any working directory.

## Input path policy

- **Existing analysis inputs** (no flags needed): `04_hotspot_analysis/data/methylated_promoter_hotspots.bed` and `04_hotspot_analysis/data/input_promoters_coordinates_methylation.tsv`.
- **Shared GSC annotation** (no flag needed): `/gsc/www/bcgsc.ca/downloads/promoter_methylation_nanopore/annotations/genome/hg38.chrom.sizes`.
- **Downloaded external datasets** (remain outside GitHub): original GSC `.../promoter_clustering/ChrmHMM_downloads`, `CTCF_binding_sites`, `H3K27me3`, `H3K4me3`, and `TAD_downloads`. These paths remain the defaults. Override their common parent with `--reference-root /path/to/promoter_clustering` if the downloads have been moved. The flag does not copy anything.
- **Results** (no flags needed): `04_hotspot_analysis/results/`, with the original subfolder names expected by the final plotting script.

The `--hotspots`, `--promoters`, `--chrom-sizes`, `--outdir`, `--permutations`, and `--seed` options remain available for non-default environments. You should **not need them on GSC** when your data files are in the existing `data/` directory and the downloaded references remain in their original locations.

## Execute on GSC

From `04_hotspot_analysis/` (or by absolute script path):

```bash
python analyze_hotspot_ChromHMM.py
python summarize_ChromHMM_across_epigenomes.py
python analyze_hotspot_CTCF.py
python analyze_hotspot_histone_marks.py
python analyze_hotspot_TAD_overlap.py
# Run the already validated CpG analysis if its output is not yet present:
python analyze_hotspot_CpG_features.py
python plot_hotspot_functional_features.py
```

The figure script reads six expected summary TSV files from `results/` and writes the figure to `results/reviewer_functional_feature_figures/`. If you already have validated CpG output, you can reuse it in `results/CpG_sequence_analysis/` rather than rerunning it.

## External datasets and dependencies

- ChromHMM: 10 Roadmap epigenomes (18-state BED.gz); NumPy, pandas.
- CTCF: 7 ENCODE peak datasets; NumPy, pandas.
- Histone: H3K27me3 and H3K4me3 bigWigs for 8 contexts; NumPy, pandas, pyBigWig.
- TAD: 7 TAD maps; NumPy, pandas.
- ChromHMM cross-epigenome summary: SciPy, Matplotlib, NumPy, pandas.
- Figure: Matplotlib, NumPy, pandas.

`reference_manifest.tsv` records accession/file identifiers from the original scripts. The source reference files **are not included** in the ZIP or repository. `.gitignore` additionally excludes common large reference formats and `results/`.

## Validation and reproducibility

The original analysis logic, permutation counts (10,000) and random seed (12345) were retained. This package has passed Python syntax compilation, but the full GSC datasets are not accessible here, so numerical results have **not** been independently validated. Compare the generated TSVs and plots with the original outputs before publication. The cCRE and older single-context analyses are omitted.
