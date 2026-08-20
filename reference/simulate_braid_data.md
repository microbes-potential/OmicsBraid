# Simulate multi-omics data with known braid patterns

Generates matched sample-level data for method development and
validation. Residuals are correlated across omic layers for each entity,
allowing the covariance, heterogeneity, equivalence, and classification
procedures to be tested against known cross-layer effects.

## Usage

``` r
simulate_braid_data(
  n_per_group = 60L,
  n_reference = NULL,
  n_comparison = NULL,
  omics = c("RNA", "Protein", "Metabolite"),
  patterns = NULL,
  rho = 0.4,
  missing_rate = 0,
  modality_missing_rate = 0,
  residual_distribution = c("normal", "t"),
  t_df = 5,
  seed = 1L
)
```

## Arguments

- n_per_group:

  Default samples per group when \`n_reference\` and \`n_comparison\`
  are not supplied.

- n_reference:

  Optional reference-group sample size.

- n_comparison:

  Optional comparison-group sample size.

- omics:

  Ordered omic names.

- patterns:

  Named list of true standardized mean shifts, one numeric vector per
  entity.

- rho:

  Scalar equicorrelation or an omic-by-omic residual correlation matrix.

- missing_rate:

  Independent value-level missingness probability.

- modality_missing_rate:

  Probability of an entire subject modality being absent; scalar or
  named by omic.

- residual_distribution:

  \`"normal"\` or heavy-tailed \`"t"\` residuals.

- t_df:

  Degrees of freedom for t residuals; must exceed 2 for finite variance.

- seed:

  Random seed.

## Value

List with \`data\` (\`omics_braid_data\`) and \`truth\` table.
