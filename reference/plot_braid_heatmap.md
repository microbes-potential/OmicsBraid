# Plot a braid heatmap across entities and omics

Plot a braid heatmap across entities and omics

## Usage

``` r
plot_braid_heatmap(
  result,
  entities = NULL,
  omic_order = result$settings$omic_order
)
```

## Arguments

- result:

  An \`omics_braid_result\`.

- entities:

  Optional entities; default top 25 by integrated adjusted p-value.

- omic_order:

  Optional omic order.

## Value

A ggplot object.
