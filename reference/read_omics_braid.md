# Read OmicsBraid input files

Read OmicsBraid input files

## Usage

``` r
read_omics_braid(
  metadata_file,
  assay_files,
  sample_id = "sample_id",
  sep = NULL,
  annotation_file = NULL
)
```

## Arguments

- metadata_file:

  CSV/TSV file containing sample metadata.

- assay_files:

  Named character vector or named list mapping omic names to CSV/TSV
  files.

- sample_id:

  Metadata sample identifier column.

- sep:

  Separator. If \`NULL\`, inferred from file extension
  (\`.tsv\`/\`.txt\` -\> tab; otherwise comma).

- annotation_file:

  Optional annotation CSV/TSV file.

## Value

An \`omics_braid_data\` object.
