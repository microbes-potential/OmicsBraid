# Test ordered cross-omic effect trajectories with generalized least squares

Fits a covariance-aware linear trajectory to standardized effects across
an explicitly ordered set of omic layers. Effects are aligned to the
dominant observed direction before fitting, so a positive slope
represents increasing absolute effect magnitude (amplification) and a
negative slope represents decreasing magnitude (attenuation). Three
practical hypotheses are evaluated: a meaningfully positive slope, a
meaningfully negative slope, and practical equivalence of the slope to a
flat trajectory within \`\[-trajectory_margin, +trajectory_margin\]\`.

## Usage

``` r
test_braid_trend(
  effects,
  covariance = NULL,
  omic_order,
  trajectory_margin = 0.15,
  alpha = 0.05,
  min_omics = 2L,
  p_adjust = "BH"
)
```

## Arguments

- effects:

  Data frame containing \`entity\`, \`omic\`, \`effect\`, and \`se\`.

- covariance:

  Optional output of \`bootstrap_effect_covariance()\` or a named list
  of entity-specific covariance matrices. If omitted, layer estimates
  are treated as independent for this calculation.

- omic_order:

  Ordered character vector describing the layer trajectory.

- trajectory_margin:

  Smallest meaningful change in standardized effect per one-layer
  transition. A scalar greater than zero.

- alpha:

  Local significance level used to define the trend state.

- min_omics:

  Minimum observed layers required.

- p_adjust:

  Multiple-testing method used for confirmatory adjusted trend p-values
  across entities. Local states remain the default for braid geometry.

## Value

One row per entity containing the aligned GLS slope, uncertainty,
practical trend tests, and local/adjusted trajectory states.

## Details

The trend test is intended for ordered layers when
attenuation/amplification is scientifically meaningful. It does not
establish causality or temporal direction.
