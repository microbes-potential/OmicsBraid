# Classify cross-omic braid patterns using inferential evidence

Braid labels are deliberately conservative. Opposite statistically
supported layer directions confirm inversion. Buffering and emergence
require practical equivalence in the appropriate downstream/upstream
layers. Concordance, attenuation, and amplification require all observed
layers to support one direction and use a covariance-aware GLS
trajectory test. If the joint-null omnibus test is not rejected and
practical equivalence is not established, the result is labelled
\`no_detectable_effect\` rather than incorrectly claiming equivalence.

## Usage

``` r
classify_braids(
  equivalence,
  omic_order,
  covariance = NULL,
  integrated = NULL,
  trend = NULL,
  trajectory_margin = 0.15,
  alpha = 0.05,
  min_slope = NULL
)
```

## Arguments

- equivalence:

  Output of \`test_equivalence()\`.

- omic_order:

  Ordered character vector describing the biological/display order of
  omics.

- covariance:

  Optional covariance object used when a trend table must be computed
  internally.

- integrated:

  Optional output of \`integrate_effects()\`. Supplying it allows
  \`no_detectable_effect\` to be distinguished from generic uncertainty.

- trend:

  Optional output of \`test_braid_trend()\`. If absent it is computed.

- trajectory_margin:

  Smallest meaningful effect change per one-layer transition for
  attenuation/amplification.

- alpha:

  Local inferential significance level.

- min_slope:

  Deprecated alias for \`trajectory_margin\` retained for early
  OmicsBraid prototypes.

## Value

One row per entity containing confirmatory and suggestive labels plus
trajectory diagnostics.

## Details

A separate \`suggestive_pattern\` is derived from effect geometry only.
It can be useful when the inferential pattern is unresolved because the
data are too imprecise, but it is not a confirmatory conclusion.
