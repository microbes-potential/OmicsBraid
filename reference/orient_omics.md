# Harmonize the sign orientation of omic-layer effects

Some omic measurements have an interpretation whose natural direction is
opposite to an activity/abundance scale used in other layers. This
helper multiplies selected omic effects by +1 or -1 and applies the
corresponding sign transformation to covariance matrices. Use it only
when the orientation is scientifically justified; it must not be used to
force apparent agreement.

## Usage

``` r
orient_omics(effects, covariance = NULL, orientation)
```

## Arguments

- effects:

  Effect table with \`omic\` and \`effect\` columns.

- covariance:

  Optional OmicsBraid covariance object or named covariance list.

- orientation:

  Named numeric vector with one value (+1 or -1) per omic to be
  re-oriented. Omics not named default to +1.

## Value

A list with oriented \`effects\` and \`covariance\`.
