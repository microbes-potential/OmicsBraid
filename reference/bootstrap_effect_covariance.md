# Estimate cross-omic effect covariance by matched-subject bootstrap

Biological subjects are resampled as whole units, preserving their
matched measurements across omics. Bootstrap correlations are estimated
for each entity and combined with analytic marginal standard errors from
\`estimate_effects()\`: V = diag(SE) R_boot diag(SE). This stabilizes
marginal uncertainty while retaining empirically estimated cross-layer
dependence.

## Usage

``` r
bootstrap_effect_covariance(
  data,
  group,
  reference,
  comparison,
  effects = NULL,
  entities = NULL,
  B = 500L,
  seed = 1L,
  min_n = 3L,
  shrinkage = 0.05,
  min_complete = 50L,
  stratified = TRUE
)
```

## Arguments

- data:

  An \`omics_braid_data\` object.

- group:

  Metadata group column.

- reference:

  Reference group.

- comparison:

  Comparison group.

- effects:

  Optional output of \`estimate_effects()\`.

- entities:

  Optional entities to bootstrap. By default, entities present in at
  least two omics.

- B:

  Number of subject-level bootstrap replicates.

- seed:

  Random seed.

- min_n:

  Minimum observations per group within an omic.

- shrinkage:

  Correlation shrinkage toward the identity in \[0,1\].

- min_complete:

  Minimum usable bootstrap pairs to estimate a correlation.

- stratified:

  Logical; resample subjects separately within reference and comparison
  groups. This preserves the observed group sizes and is recommended for
  fixed two-group designs.

## Value

Object of class \`omics_braid_covariance\` containing a covariance
matrix per entity.
