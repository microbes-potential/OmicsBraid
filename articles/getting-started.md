# Getting started with OmicsBraid

## Purpose

OmicsBraid analyzes the **same biological contrast** across
scientifically ordered molecular layers while retaining layer-specific
uncertainty and matched-subject dependence.

## Minimal known-truth example

``` r

library(OmicsBraid)
sim <- simulate_braid_data(
  n_per_group = 40,
  omics = c("RNA", "Protein", "Phosphoprotein"),
  rho = 0.4,
  seed = 42
)
fit <- run_omics_braid(
  sim$data,
  group = "group",
  reference = "Control",
  comparison = "Disease",
  omic_order = c("RNA", "Protein", "Phosphoprotein"),
  bootstrap_B = 300,
  equivalence_margin = 0.30,
  trajectory_margin = 0.15,
  pattern_draws = 1000,
  seed = 42
)
```

## Inspect results

``` r

fit$effects
fit$integrated
fit$equivalence
fit$trend
fit$classification
braid_results_table(fit)
```

`fit$integrated` contains the GLS consensus, omnibus evidence, and
cross-omic heterogeneity. The classification layer should be interpreted
together with `pattern_status` and the underlying effect vector rather
than as an isolated label.

## Visualize one entity

``` r

plot_evidence_forest(fit, "inversion")
plot_effect_braid(fit, "inversion")
```

For a genome-wide overview:

``` r

plot_concordance_map(fit)
plot_braid_heatmap(fit)
```

## Real data

Create an `omics_braid_data` object from analysis-ready assay matrices
and metadata. Column names must identify samples and metadata must
contain the same sample IDs.

``` r

obj <- omics_braid_data(
  assays = list(RNA = rna, Protein = protein, Phosphoprotein = phospho),
  metadata = metadata,
  sample_id = "sample_id"
)
validate_omics_braid_data(obj)
```

Do not feed raw sequencing counts or raw mass-spectrometry files
directly into OmicsBraid.
