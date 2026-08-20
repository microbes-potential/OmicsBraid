# Empirically calibrate OmicsBraid omnibus and heterogeneity tests

Provides resampling-based p-values for the two cross-omic quadratic
tests. The omnibus test can be calibrated by matched-subject label
permutation or by a centered matched-subject bootstrap. The
heterogeneity test can be calibrated by a raw-data null-shift matched
bootstrap (recommended) or by an effect-level centered bootstrap under
the fitted common-effect null.

## Usage

``` r
empirical_omics_tests(
  data,
  group,
  reference,
  comparison,
  effects = NULL,
  entities = NULL,
  B = 499L,
  seed = 1L,
  min_n = 3L,
  min_omics = 2L,
  min_complete = 100L,
  omnibus_method = c("permutation", "centered_bootstrap"),
  heterogeneity_method = c("null_shift_bootstrap", "centered_bootstrap"),
  orientation = NULL,
  p_adjust = "BH"
)
```

## Arguments

- data:

  An \`omics_braid_data\` object containing sample-level assays.

- group:

  Metadata column containing the two groups.

- reference:

  Reference-group label.

- comparison:

  Comparison-group label.

- effects:

  Optional layer-specific effect table. If supplied after
  \`orient_omics()\`, pass the same \`orientation\` so resampled effects
  receive the identical sign transformation.

- entities:

  Optional entities to calibrate. By default, entities observed in at
  least \`min_omics\` layers are used.

- B:

  Number of resampling replicates for each empirical null.

- seed:

  Random seed.

- min_n:

  Minimum observations per group within an omic.

- min_omics:

  Minimum omic layers per entity.

- min_complete:

  Minimum complete resampling draws required for a p-value.

- omnibus_method:

  Either \`"permutation"\` or \`"centered_bootstrap"\`. Permutation is
  appropriate for the global null in an exchangeable two-group design.
  Centered bootstrap is a nonparametric alternative.

- heterogeneity_method:

  \`"null_shift_bootstrap"\` (recommended robust calibration) or
  \`"centered_bootstrap"\`. Ordinary label permutation is not used
  because the heterogeneity null permits a common non-zero effect.

- orientation:

  Optional named +1/-1 vector applied to the resampled layer effects.
  This must match any scientific orientation already applied to
  \`effects\`.

- p_adjust:

  Multiple-testing method for empirical p-values.

## Value

A data frame containing asymptotic-independent empirical omnibus and
heterogeneity p-values, empirical critical values, and resampling
diagnostics.

## Details

The resampling is performed at the biological-subject level: all
available omic measurements belonging to a subject remain linked. This
preserves the cross-omic dependence that would be destroyed by shuffling
individual assay matrices independently.

Empirical calibration is intended as a robust alternative when the
chi-square reference distributions used by \`integrate_effects()\` may
be inaccurate, for example under heavy-tailed sampling distributions.
The asymptotic statistics remain available and are not overwritten by
this function.
