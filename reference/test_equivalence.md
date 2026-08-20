# Test practical equivalence to a negligible effect region

Performs two one-sided tests (TOST) for whether each standardized effect
lies within \`\[-margin, +margin\]\`. It reports both local
(entity-specific) and multiplicity-adjusted inferential states. Local
states are recommended for describing the geometry of an individual
braid because they do not change merely when unrelated entities are
added to the analysis; adjusted states are retained for confirmatory
screening across many entities.

## Usage

``` r
test_equivalence(
  effects,
  margin = 0.3,
  alpha = 0.05,
  p_adjust = "BH",
  state_basis = c("local", "adjusted")
)
```

## Arguments

- effects:

  Effect table from \`estimate_effects()\` or compatible summary
  statistics.

- margin:

  Smallest effect size of interest on the standardized-effect scale;
  scalar or named by omic.

- alpha:

  Significance level.

- p_adjust:

  Multiple-testing adjustment applied separately within each omic.

- state_basis:

  Which inferential state is copied to the legacy \`state\` column:
  \`"local"\` (default) or \`"adjusted"\`.

## Value

Effect table with TOST p-values, local/adjusted difference p-values, and
both local and multiplicity-adjusted practical states.
