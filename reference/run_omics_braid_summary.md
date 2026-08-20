# Run OmicsBraid from externally estimated summary statistics

Run OmicsBraid from externally estimated summary statistics

## Usage

``` r
run_omics_braid_summary(
  effects,
  covariance = NULL,
  omic_order,
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

- effects:

  Data frame with at least \`entity\`, \`omic\`, \`effect\`, and \`se\`.

- covariance:

  Optional named list of entity-specific covariance matrices.

- omic_order:

  Ordered omics.

- orientation:

  Optional named +1/-1 vector to harmonize omic effect directions.

- equivalence_margin:

  Smallest effect size of interest.

- alpha:

  Significance level.

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

  Monte Carlo draws.

- seed:

  Random seed.

## Value

\`omics_braid_result\` without sample-level data.
