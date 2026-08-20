# OmicsBraid

**Covariance-aware inference of cross-omic effect trajectories in matched multi-omics studies**

OmicsBraid is an R package for estimating, testing, classifying, and visualizing how a biological contrast changes across ordered molecular layers. Rather than reducing multi-omics results to a single pooled score, OmicsBraid keeps the **magnitude, direction, covariance, heterogeneity, practical equivalence, and ordered trajectory** of layer-specific effects explicit.

**Frozen manuscript release:** `v0.2.2`

> Statistical core status: **frozen for manuscript use.** The `R/` source in this repository is byte-for-byte identical to the validated OmicsBraid v0.2.2 core used in the simulation, real-background semi-synthetic, comparator/ablation, CPTAC-GBM, and independent CPTAC-LUAD analyses.

## What OmicsBraid answers

For an ordered effect vector such as

```text
RNA  ->  Protein  ->  Phosphoprotein
```

OmicsBraid separates four questions:

1. **Layer-specific effects:** what is the standardized effect in each omic?
2. **Global evidence:** is the joint multi-omic effect vector different from zero?
3. **Consensus:** is there a common covariance-aware effect across layers?
4. **Trajectory/heterogeneity:** does the effect attenuate, amplify, buffer, emerge, invert, or remain concordant across the ordered layers?

The framework is not restricted to RNA, protein, and phosphoprotein data. Any scientifically defensible ordered set of analysis-ready molecular layers can be used.

## Statistical framework

OmicsBraid implements:

- Hedges' g layer-specific standardized effects;
- matched-subject bootstrap estimation of cross-omic dependence;
- generalized least-squares (GLS) consensus effects;
- covariance-aware multivariate omnibus inference;
- `Q_omics` cross-omic heterogeneity inference;
- TOST practical-equivalence testing;
- covariance-aware ordered GLS trajectory testing;
- hierarchical braid classification with confirmatory and suggestive states;
- frequentist Monte Carlo propagation of braid geometry uncertainty;
- subject-level permutation calibration of the omnibus test;
- null-shift matched-subject bootstrap calibration of `Q_omics`;
- analytic confidence intervals with BCa sensitivity support;
- Omics Evidence Forest, Effect Braid, braid heatmap, and concordance-map visualizations.

A near-zero GLS consensus is **not** interpreted as absence of biology when the multivariate effect or heterogeneity is strong. This is particularly important for cross-layer inversions, where signed pooling can cancel opposing effects.

## Installation

### Development/manuscript release from GitHub

The public release metadata is already configured for the `microbes-potential/OmicsBraid` repository. Users can install the frozen manuscript release with:

```r
install.packages("remotes")
remotes::install_github("microbes-potential/OmicsBraid", ref = "v0.2.2")
```

For a local source checkout:

```r
install.packages(".", repos = NULL, type = "source")
```

## Quick start

The package includes a known-truth simulator so users can try the complete workflow without downloading external data.

```r
library(OmicsBraid)

sim <- simulate_braid_data(
  n_per_group = 40,
  omics = c("RNA", "Protein", "Phosphoprotein"),
  rho = 0.4,
  seed = 42
)

fit <- run_omics_braid(
  data = sim$data,
  group = "group",
  reference = "Control",
  comparison = "Disease",
  omic_order = c("RNA", "Protein", "Phosphoprotein"),
  bootstrap_B = 300,
  bootstrap_shrinkage = 0.05,
  empirical_tests = FALSE,
  ci_method = "analytic",
  integrated_ci_method = "analytic",
  equivalence_margin = 0.30,
  trajectory_margin = 0.15,
  pattern_draws = 1000,
  seed = 42
)

braid_results_table(fit)
plot_evidence_forest(fit, "inversion")
plot_effect_braid(fit, "inversion")
```

For targeted robust empirical calibration:

```r
emp <- empirical_omics_tests(
  data = sim$data,
  group = "group",
  reference = "Control",
  comparison = "Disease",
  effects = fit$effects,
  B = 999,
  seed = 42,
  omnibus_method = "permutation",
  heterogeneity_method = "null_shift_bootstrap",
  p_adjust = "BH"
)
```

## Primary braid classes

| Class | Interpretation |
|---|---|
| `concordant_increase` | supported increase across ordered layers |
| `concordant_decrease` | supported decrease across ordered layers |
| `attenuation` | same-direction effect becomes meaningfully weaker across layers |
| `amplification` | same-direction effect becomes meaningfully stronger across layers |
| `buffering` | upstream effect with downstream practical equivalence |
| `emergence` | upstream practical equivalence with a downstream effect |
| `inversion` | supported change in effect direction across layers |
| `null_equivalent` | effects are supported as practically negligible |
| `no_detectable_effect` | insufficient evidence for a global effect; not equivalent to practical equivalence |
| `uncertain` | evidence does not support a stable confirmatory trajectory label |
| `insufficient` | too little usable cross-omic information |

OmicsBraid deliberately distinguishes **confirmed/direction-confirmed** patterns from **suggestive/unresolved** geometry. Failure to reject a null hypothesis is never treated as evidence of equivalence.

## Input

Each assay is a numeric matrix with rows as entities and columns as biological samples. Sample IDs link assays to metadata.

```r
obj <- omics_braid_data(
  assays = list(
    RNA = rna_matrix,
    Protein = protein_matrix,
    Phosphoprotein = phosphoprotein_matrix
  ),
  metadata = metadata,
  sample_id = "sample_id"
)
```

The package expects **analysis-ready** omic values. Raw RNA-seq counts and raw mass-spectrometry files require assay-specific preprocessing before OmicsBraid.

Small CSV/TSV templates are distributed in `inst/extdata/` and are installed with the package.

## Recommended inference

For routine reporting, OmicsBraid v0.2.2 retains analytic layer and consensus confidence intervals. When non-normality is a concern, BCa layer intervals can be used as a sensitivity analysis. For robust hypothesis inference, the validated empirical procedures are:

- **global omnibus:** subject-level group-label permutation;
- **cross-omic heterogeneity:** null-shift matched-subject bootstrap.

Finite empirical p-values are not automatically suitable for genome-wide BH-FDR when the number of resamples is too small to provide adequate p-value resolution. The manuscript applications therefore used continuous asymptotic p-values for genome-wide screening and empirical tests as pre-declared secondary robustness checks on candidate subsets.

## Scope and safeguards

OmicsBraid v0.2.2 is intended for research use with analysis-ready bulk/sample-level multi-omics data and two independent biological groups. It supports matched subjects across layers and partial modality missingness.

The current release does **not** silently approximate paired/repeated-measures designs, survival outcomes, continuous exposures, more than two comparison groups, or causal molecular propagation. Omic ordering, orientation, equivalence margins, and trajectory margins must be scientifically justified before outcome-driven interpretation.

## Documentation

- `vignette("OmicsBraid-introduction", package = "OmicsBraid")`
- `vignette("getting-started", package = "OmicsBraid")`
- `vignette("statistical-framework", package = "OmicsBraid")`
- `vignette("braid-classification", package = "OmicsBraid")`
- `vignette("robust-inference", package = "OmicsBraid")`
- `vignette("complete-workflow", package = "OmicsBraid")`

The pkgdown website is configured in `_pkgdown.yml`; GitHub Actions can publish it automatically after the repository is created.

## Validation and reproducibility

The statistical core used here was validated through:

- broad known-truth simulation;
- robust Type-I calibration under normal and heavy-tailed settings;
- confidence-interval sensitivity experiments;
- CPTAC real-background semi-synthetic validation preserving observed covariance and missingness;
- comparator/ablation experiments;
- a genuine CPTAC-GBM application; and
- an independent CPTAC-LUAD application with ranked pathway validation.

The manuscript-scale analysis scripts and derived result tables belong in the separate **`OmicsBraid-paper`** reproducibility repository, not in this software repository. Raw CPTAC data are not redistributed.

## Citation

```r
citation("OmicsBraid")
```

The repository also contains `CITATION.cff` for GitHub/Zenodo metadata. Update the preferred manuscript citation after the article receives a DOI.

## Versioning

`v0.2.2` is the frozen manuscript analysis release. Future software improvements should use a new version and must not overwrite or silently alter the `v0.2.2` tag.

## License

MIT License. See `LICENSE` and `LICENSE.md`.

## Disclaimer

OmicsBraid is research software and is not intended for clinical decision-making.
