# Plot an Effect Braid

Plot an Effect Braid

## Usage

``` r
plot_effect_braid(
  result,
  entity,
  omic_order = result$settings$omic_order,
  show_ci = TRUE,
  show_equivalence_region = TRUE
)
```

## Arguments

- result:

  An \`omics_braid_result\`.

- entity:

  Entity/pathway to display.

- omic_order:

  Optional omic order.

- show_ci:

  Show confidence intervals.

- show_equivalence_region:

  Shade the negligible-effect region.

## Value

A ggplot object.
