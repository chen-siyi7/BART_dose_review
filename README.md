# BART_dose_review

An R package for **Checking and Explaining BART-Guided Dose Recommendations**, by Siyi Chen, School of Public Health, LSU Health Sciences Center.

The package checks which dose recommendations are supported by numerical intervals, finds shorter reports that support the same recommendation, and examines changes that could permit another dose.

The GitHub repository is `BART_dose_review`; the R package name is `bartDoseReview`.

## Install

Use R 4.4 or later. From the directory containing the downloaded repository:

```r
install.packages(c("digest", "gmp", "jsonlite"))
install.packages("BART_dose_review", repos = NULL, type = "source")
library(bartDoseReview)
```

Clinical replay also uses `posterior` and `readxl`; figures use `ggplot2`. Fresh BART fits use version 2.9.9 of `BART`.

## Example

```r
library(bartDoseReview)
gates <- rep(list(rbind(c(0.7, 0.9), c(0.1, 0.2))), 3)
scores <- rbind(c(0.1, 0.2), c(0.7, 0.8), c(0.4, 0.5))
interval_actions("Final", gates, scores)
```

The result is `"dose:2"`. Use `help(package = "bartDoseReview")` for function documentation.

## Reproduce the analyses

Use the separate full reproducibility bundle, which contains the data, saved draws and run records. Set `bundle` to its directory containing `MANIFEST.json`, `replay_all.R` and `data/`, and choose a new output directory:

```r
bundle <- "/path/to/reproducibility"
replay_results(bundle, "/path/to/new/replay_output")
```

This verifies the bundle and its source version before reconstructing the saved results. Fresh analyses and figures can be run with:

```r
run_analysis("feasibility", bundle, "/path/to/new/feasibility")
run_analysis("figures", bundle, "/path/to/new/figures")
```

`?run_analysis` lists the supported studies. The clinical source data are from [Waldbaum et al. (2020)](https://doi.org/10.1371/journal.pone.0232394); original files and attribution remain in the separate bundle.

## Layout

- `R/`: package functions.
- `man/`: function help files.
- `inst/scripts/`: command-line entry point.
- `inst/studies/` and `inst/clinical/`: the study and clinical source scripts used by the full bundle. Run these workflows through `run_analysis()` so their input paths are resolved in the bundle.
- `tests/`: checks for decisions, minimum reports and reference calculations.
- `tools/verify.R`: source and package-layout checks.

All R code is comment-free. Manuscripts, datasets, saved results, figures and historical source copies are kept outside the package. `DESCRIPTION` and `NAMESPACE` provide the package metadata and exports.

Repository: [chen-siyi7/BART_dose_review](https://github.com/chen-siyi7/BART_dose_review).
