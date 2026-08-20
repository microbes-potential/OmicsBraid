# Bootstrap confidence intervals for GLS consensus effects

Uses matched subject-bootstrap layer-effect draws from
\`bootstrap_effect_covariance()\` and the fitted covariance structure to
form a bootstrap distribution of the GLS consensus effect. The
hypothesis tests and analytic standard errors remain unchanged; these
intervals are intended for robust uncertainty reporting and coverage
validation.

## Usage

``` r
bootstrap_consensus_intervals(
  effects,
  bootstrap,
  integrated = NULL,
  method = c("percentile", "basic"),
  conf_level = 0.95,
  min_boot = 100L,
  min_omics = 2L
)
```

## Arguments

- effects:

  Layer-specific effect table.

- bootstrap:

  An \`omics_braid_covariance\` object.

- integrated:

  Optional output of \`integrate_effects()\`.

- method:

  \`"percentile"\` or \`"basic"\`.

- conf_level:

  Confidence level.

- min_boot:

  Minimum usable bootstrap consensus draws.

- min_omics:

  Minimum finite omic layers required in an individual bootstrap
  replicate.

## Value

One row per entity with bootstrap consensus interval diagnostics.
