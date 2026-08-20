# Plot concordance versus integrated significance

Plot concordance versus integrated significance

## Usage

``` r
plot_concordance_map(
  result,
  label_top = 0L,
  metric = c("evidence_direction_agreement", "direction_agreement", "i2_consistency"),
  significance = c("omnibus", "consensus")
)
```

## Arguments

- result:

  An \`omics_braid_result\`.

- label_top:

  Number of most significant entities to label.

- metric:

  Concordance metric: evidence-qualified directional agreement
  (default), raw weighted directional agreement, or descriptive I2-based
  consistency.

- significance:

  Evidence axis: multivariate omnibus (default) or GLS consensus-effect
  significance.

## Value

A ggplot object.
