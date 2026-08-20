#' OmicsBraid: cross-omics effect concordance and heterogeneity
#'
#' OmicsBraid is a research-oriented framework for cross-layer inference from
#' matched bulk multi-omics studies. Its core estimand is the vector of
#' standardized group effects for the same biological entity/pathway across
#' ordered omics layers. It is intentionally not a raw-data preprocessing or
#' latent-factor integration package. Robust empirical calibration is available
#' for the omnibus and heterogeneity tests when asymptotic chi-square reference
#' distributions are questionable.
#'
#' @keywords internal
"_PACKAGE"

# Declare ggplot2 non-standard evaluation symbols for R CMD check.
utils::globalVariables(c(
  "concordance", "conf_high", "conf_low", "effect", "entity", "label",
  "minus_log10_p", "n_omics", "omic", "state", "type",
  "xmax", "xmin", "ymax", "ymin"
))
