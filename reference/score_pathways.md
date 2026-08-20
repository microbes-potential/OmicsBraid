# Score pathways within each omic layer

This convenience function z-standardizes each feature across samples
within an omic and then aggregates features assigned to the same
pathway. OmicsBraid's inferential core can also accept externally
computed pathway/activity scores, which is recommended when a
domain-specific scoring method is preferred.

## Usage

``` r
score_pathways(
  data,
  mapping,
  method = c("mean_z", "median_z"),
  min_features = 3L,
  center = TRUE,
  scale = TRUE
)
```

## Arguments

- data:

  An \`omics_braid_data\` object.

- mapping:

  Data frame with columns \`omic\`, \`feature_id\`, and \`pathway\`.

- method:

  Aggregation method: \`"mean_z"\` or \`"median_z"\`.

- min_features:

  Minimum mapped features per pathway within an omic.

- center:

  Logical; center feature values before aggregation.

- scale:

  Logical; scale feature values before aggregation.

## Value

An \`omics_braid_data\` object whose assay rows are pathways.
