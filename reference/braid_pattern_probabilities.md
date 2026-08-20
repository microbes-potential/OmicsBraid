# Quantify uncertainty in effect-braid geometry

Draws from a multivariate normal approximation to the estimated
cross-omic effects and reports how often each practical geometric
pattern occurs. These frequencies propagate estimation uncertainty; they
are not Bayesian posterior probabilities and they do not replace the
confirmatory inferential braid label.

## Usage

``` r
braid_pattern_probabilities(
  effects,
  covariance = NULL,
  omic_order,
  margin = 0.3,
  n_draws = 2000L,
  trajectory_margin = 0.15,
  seed = 1L,
  min_slope = NULL
)
```

## Arguments

- effects:

  Effect table.

- covariance:

  Optional bootstrap covariance object/list. If absent, omics are
  treated as independent.

- omic_order:

  Ordered omics.

- margin:

  Equivalence/negligible-effect margin; scalar or named by omic.

- n_draws:

  Number of Monte Carlo draws per entity.

- trajectory_margin:

  Practical trajectory threshold per layer transition.

- seed:

  Random seed.

- min_slope:

  Deprecated alias for \`trajectory_margin\`.

## Value

Long data frame of geometric pattern frequencies per entity.
