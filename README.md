# Triangulating evidence to detect signatures of stabilizing selection acting on protein levels in humans

Code - manuscript by **Mihaela-Diana Zanoaga, Théo Cavinato, Thibault Latrille, Tabea Schoeler, Orsolya Liska, Balázs Papp and Zoltán Kutalik** (University of Lausanne, SIB, Unisanté, Biological Research Centre Szeged).

## Overview

We investigate stabilizing selection on the circulating levels of 2'923 plasma proteins (UK Biobank Pharma Proteomics Project, Olink) in up to 337'386 UK Biobank participants by combining two complementary signatures of selection:

1. **Frequency-dependent genetic architecture**: variants with larger effects on protein levels are expected to be rarer under stabilizing selection.
   - **LDpred2-auto** on common variants (trans-pQTL summary statistics, cis region ±1 Mb excluded) → selection parameter α
   - **Rare loss-of-function burden** (MAF < 10⁻³) with a maximum-likelihood model → α_b
2. **Contemporary fitness associations**: quadratic regressions of three fitness proxies (number of children, a composite health score, parental longevity) on protein levels → quadratic coefficient *c*, where *c* < 0 with an optimum within ±1 SD indicates stabilizing selection.

The resulting estimates are then compared with each other and with independent measures of evolutionary and functional constraint (SNP-heritability, pLI, LOEUF, cross-species transcript evolutionary rate, LASSO R² on complex traits, dN/dS, ω_A, Gene Damage Index, Moran's I).

## Repository structure

```
Preprocess/           sample QC, covariates, protein levels and fitness proxies
Regression/           quadratic regressions of fitness proxies on protein levels
LDpred2/              LDpred2-auto on common variants
Burden/               likelihood model on rare loss-of-function burden estimates
Constraint_metrics/   independent constraint measures 
Triangulation/        integration of all selection estimates: patterns, concordance and overlap
```

Scripts are numbered in the order in which they are run. Steps in `Preprocess/`, `Regression/`, `LDpred2/` (01–04), `Burden/` and the LASSO script use individual-level UK Biobank data and were run on the Urblauna HPC cluster (University of Lausanne). All other steps run locally on the resulting summary tables.

### `Preprocess/`

| Script | Description |
|---|---|
| `01_sample_qc.R` | Unrelated, White British participants; redacted and withdrawn samples removed |
| `02_covariates_complex_traits.R` | Covariates (age, age², sex, genotyping batch, PC1–40) and 48 complex traits |
| `03_protein_levels.R` | Olink protein levels at baseline for the QC-passed sample (n = 36,185) |
| `04_health_score_step1.R`, `05_health_score_step2.R` | Composite health score |
| `06_parental_longevity.R` | Parental longevity |
| `07_educational_attainment.R` | Educational attainment |
| `08_number_of_children.R` | Number of children|
| `09_tea_intake_negative_control.R` | Tea intake, used as a negative-control phenotype |

### `Regression/`

| Script | Description |
|---|---|
| `01_quadratic_regression_fitness_proxies.R` | Covariate-adjusted fitness proxy ~ protein + protein² for every protein and proxy |
| `02_quadratic_regression_tea_intake.R` | Same model with tea intake as outcome (negative control) |
| `03_process_regression_results.R` | FDR, inferred optimum (−b / 2c) and flags for optimum within ±0.5 / ±1 SD |

### `LDpred2/`

| Script | Description |
|---|---|
| `01_ldpred2_auto_trans.R`, `02_combine_ldpred2_trans.R` | LDpred2-auto per protein on HapMap3+ variants (MAF > 0.01), cis region excluded: main analysis |
| `03_ldpred2_auto_cis_trans.R`, `04_combine_ldpred2_cis_trans.R` | Same, keeping cis variants (sensitivity analysis) |
| `05_process_ldpred2_cis_trans.R`, `06_process_ldpred2_trans.R` | Process LDpred2 results |
| `07_compare_cis_trans.R` | Comparison of α, h² and polygenicity estimated with and without cis variants |

### `Burden/`

| Script | Description |
|---|---|
| `01_prepare_burden_per_protein.R` | Per-protein gene-burden estimates (high-confidence LoF, MAF < 10⁻³) |
| `02_burden_likelihood_alpha.R` | Maximum-likelihood estimation of α_b |

### `Constraint_metrics/`

| Script | Description |
|---|---|
| `lasso_r2_complex_traits.R` | Cross-validated R² of a LASSO regression of each protein on the complex traits |
| `dnds_human_mouse.R` | Human–mouse dN/dS  |
| `transcript_evolutionary_rate.R` | Meta-analysis of tissue-specific transcript evolutionary rates (Liu et al. 2023) |

pLI and LOEUF (gnomAD v4.1), ω_A (Latrille et al. 2023), Gene Damage Index (Itan et al. 2015) and Moran's I were dowloaded or obtained from collaborators.

### `Triangulation/`

| Script | Description |
|---|---|
| `1.selection_categories/quadratic_selection_categories.R` | Classification of protein–fitness relationships into stabilizing, disruptive, directional non-linear and directional linear patterns |
| `1.selection_categories/complex_traits_overview.R` | Overview of the 48 complex traits by biological category |
| `2.direction_of_selection/stabilizing_vs_disruptive.R` | Number and proportion of significant negative (stabilizing) vs positive (disruptive) selection estimates per method, at nominal and FDR thresholds, including the education and tea-intake sensitivity analyses |
| `3.concordance/concordance_estimates_constraint.R` | Spearman correlations among selection estimates and with constraint metrics, for nominally and FDR-significant estimates, with the number of overlapping proteins per pair |
| `3.concordance/concordance_estimates_constraint_no_pvalue_filter.R` | Same correlations without significance filtering of the selection estimates |
| `3.concordance/significant_correlations_tables.R` | Tables of the significant correlations under the three filtering schemes |
| `4.overlap/overlap_stabilizing_proteins.R` | Overlap of proteins under stabilizing selection across the five methods |

## Filters used across analyses

- **LDpred2:** h² > 0.1, p > 0.00025, and h² and p significantly different from 0 (p < 0.05).
- **Quadratic regressions:** inferred optimum within ±1 SD of the protein distribution and h² > 0.1.
- **Burden:** proteins with estimated σ²_G = 0 (α_b not identifiable) are excluded.
- Significance is assessed at nominal p < 0.05 and at FDR < 0.05 (Benjamini–Hochberg).


## Requirements

R (≥ 4.2) with:

- **Analysis:** `bigsnpr`, `glmnet`, `orthologr`, `biomartr`, `Biostrings`
- **Data handling:** `data.table`, `dplyr`, `tidyr`, `tibble`, `readr`, `readxl`, `writexl`, `arrow`
- **Plotting:** `ggplot2`, `patchwork`, `gridExtra`, `ggrepel`, `UpSetR`

## Contact

Mihaela-Diana Zanoaga, Department of Computational Biology, University of Lausanne.
email: mihaeladiana.zanoaga@unil.ch

