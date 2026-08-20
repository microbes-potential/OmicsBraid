# Estimate layer-specific standardized effects

Estimates Hedges' g for a two-group contrast in each feature/pathway and
omic. Inputs should already be quality-controlled and normalized
appropriately for their assay technology. Hedges' g is scale-free but
does not repair poor raw preprocessing, severe censoring, or
inappropriate transformations.

## Usage

``` r
estimate_effects(
  data,
  group,
  reference,
  comparison,
  entities = NULL,
  min_n = 3L,
  conf_level = 0.95,
  p_adjust = "BH"
)
```

## Arguments

- data:

  An \`omics_braid_data\` object.

- group:

  Metadata column containing the two groups.

- reference:

  Reference group label.

- comparison:

  Comparison group label. Positive effects mean comparison \> reference.

- entities:

  Optional character vector limiting features/pathways.

- min_n:

  Minimum non-missing observations per group and omic.

- conf_level:

  Confidence level.

- p_adjust:

  Multiple-testing method applied separately within each omic.

## Value

Data frame of effect estimates and uncertainty.
