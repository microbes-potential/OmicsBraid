# Run the complete OmicsBraid workflow on sample-level data

Run the complete OmicsBraid workflow on sample-level data

## Usage

``` r
run_omics_braid(
  data,
  group,
  reference,
  comparison,
  omic_order = names(data$assays),
  entities = NULL,
  pathway_mapping = NULL,
  pathway_method = "mean_z",
  min_pathway_features = 3L,
  bootstrap_B = 500L,
  bootstrap_shrinkage = 0.05,
  empirical_tests = FALSE,
  empirical_B = 499L,
  empirical_omnibus_method = c("permutation", "centered_bootstrap"),
  empirical_heterogeneity_method = c("null_shift_bootstrap", "centered_bootstrap"),
  empirical_use_as_primary = FALSE,
  ci_method = c("analytic", "percentile", "basic", "bca"),
  integrated_ci_method = c("analytic", "percentile", "basic"),
  ci_conf_level = 0.95,
  ci_min_boot = 100L,
  orientation = NULL,
  equivalence_margin = 0.3,
  alpha = 0.05,
  state_basis = c("local", "adjusted"),
  equivalence_p_adjust = "BH",
  trajectory_margin = 0.15,
  trend_p_adjust = "BH",
  min_slope = NULL,
  pattern_draws = 2000L,
  seed = 1L
)
```

## Arguments

- data:

  An \`omics_braid_data\` object.

- group:

  Metadata group column.

- reference:

  Reference group.

- comparison:

  Comparison group.

- omic_order:

  Biological/display order of omics. Defaults to assay order.

- entities:

  Optional entities to analyze.

- pathway_mapping:

  Optional mapping passed to \`score_pathways()\`; if supplied, analysis
  is performed on pathway scores.

- pathway_method:

  Pathway scoring method.

- min_pathway_features:

  Minimum pathway features per omic.

- bootstrap_B:

  Bootstrap replicates for cross-omic covariance. Set to 0 to assume
  independence.

- bootstrap_shrinkage:

  Correlation shrinkage.

- empirical_tests:

  Logical; if \`TRUE\`, calculate additional resampling- calibrated
  omnibus and heterogeneity p-values using \`empirical_omics_tests()\`.

- empirical_B:

  Number of empirical resampling replicates.

- empirical_omnibus_method:

  \`"permutation"\` or \`"centered_bootstrap"\`.

- empirical_heterogeneity_method:

  \`"null_shift_bootstrap"\` (recommended) or \`"centered_bootstrap"\`.

- empirical_use_as_primary:

  Logical; if \`TRUE\`, empirical p-values replace asymptotic p-values
  when available for downstream omnibus evidence. The original
  asymptotic p-values are retained in separate columns.

- ci_method:

  Layer-effect confidence interval method: \`"analytic"\`,
  \`"percentile"\`, \`"basic"\`, or \`"bca"\`. Bootstrap intervals
  change interval reporting only; analytic SEs and p-values remain
  available.

- integrated_ci_method:

  Consensus-effect confidence interval method: \`"analytic"\`,
  \`"percentile"\`, or \`"basic"\`.

- ci_conf_level:

  Confidence level for analytic/bootstrap intervals.

- ci_min_boot:

  Minimum finite bootstrap draws required before a bootstrap interval
  replaces the analytic interval.

- orientation:

  Optional named +1/-1 vector to harmonize omic effect directions.

- equivalence_margin:

  Smallest effect size of interest for practical equivalence.

- alpha:

  Significance level for local difference/equivalence and braid
  inference.

- state_basis:

  Use local or multiplicity-adjusted inferential states for
  deterministic braid labels.

- equivalence_p_adjust:

  Multiple-testing method retained for adjusted equivalence/difference
  evidence.

- trajectory_margin:

  Smallest meaningful standardized-effect change per one-layer
  transition.

- trend_p_adjust:

  Multiple-testing method retained for adjusted trajectory evidence.

- min_slope:

  Deprecated alias for \`trajectory_margin\`.

- pattern_draws:

  Monte Carlo draws for geometric uncertainty propagation.

- seed:

  Random seed.

## Value

An object of class \`omics_braid_result\`.
