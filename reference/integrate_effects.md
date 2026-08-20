# Integrate effects across omic layers and quantify heterogeneity

Uses generalized least squares (GLS) to estimate a common cross-omic
effect. A generalized Cochran Q statistic tests whether the
layer-specific effects are compatible with a common effect after
accounting for their sampling covariance. A separate multivariate
Wald-type omnibus statistic tests the joint null that all layer effects
are zero; unlike the consensus effect, this test does not cancel equally
strong effects occurring in opposite directions. The reported I2-like
statistic is descriptive and should not be interpreted as literal
between-study heterogeneity because omics layers are not studies.

## Usage

``` r
integrate_effects(
  effects,
  covariance = NULL,
  min_omics = 2L,
  p_adjust = "BH",
  conf_level = 0.95
)
```

## Arguments

- effects:

  Data frame containing \`entity\`, \`omic\`, \`effect\`, and \`se\`.

- covariance:

  Optional output of \`bootstrap_effect_covariance()\` or a named list
  of covariance matrices.

- min_omics:

  Minimum omics per entity.

- p_adjust:

  Multiple-testing method for integrated and heterogeneity p-values.

- conf_level:

  Confidence level for the analytic GLS consensus interval.

## Value

Data frame of integrated effects and heterogeneity diagnostics.
