# Harmonize assay-specific feature identifiers to common entities

Maps feature IDs within each omic to a shared entity identifier for
entity-level cross-omic analysis (for example Ensembl RNA identifiers
and UniProt proteins to a common gene symbol). Many-to-one mappings are
rejected by default because collapsing isoforms/probes changes the
scientific estimand.

## Usage

``` r
harmonize_entities(
  data,
  mapping,
  collapse = c("error", "mean", "median"),
  min_mapped = 1L
)
```

## Arguments

- data:

  An \`omics_braid_data\` object.

- mapping:

  Data frame with columns \`omic\`, \`feature_id\`, and \`entity\`.

- collapse:

  How to handle multiple assay features mapping to one entity:
  \`"error"\` (default), \`"mean"\`, or \`"median"\`.

- min_mapped:

  Minimum mapped entities required per retained omic.

## Value

An \`omics_braid_data\` object with harmonized entity row names and a
mapping report stored in \`attr(x, "entity_harmonization")\`.
