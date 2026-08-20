# Construct an OmicsBraid data object

Construct an OmicsBraid data object

## Usage

``` r
omics_braid_data(assays, metadata, sample_id = "sample_id", annotation = NULL)
```

## Arguments

- assays:

  Named list of numeric matrices. Features are rows and samples are
  columns.

- metadata:

  Data frame containing one row per biological sample.

- sample_id:

  Column in \`metadata\` holding sample identifiers.

- annotation:

  Optional feature annotation data frame.

## Value

An object of class \`omics_braid_data\`.
