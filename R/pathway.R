#' Score pathways within each omic layer
#'
#' This convenience function z-standardizes each feature across samples within an
#' omic and then aggregates features assigned to the same pathway. OmicsBraid's
#' inferential core can also accept externally computed pathway/activity scores,
#' which is recommended when a domain-specific scoring method is preferred.
#'
#' @param data An `omics_braid_data` object.
#' @param mapping Data frame with columns `omic`, `feature_id`, and `pathway`.
#' @param method Aggregation method: `"mean_z"` or `"median_z"`.
#' @param min_features Minimum mapped features per pathway within an omic.
#' @param center Logical; center feature values before aggregation.
#' @param scale Logical; scale feature values before aggregation.
#' @return An `omics_braid_data` object whose assay rows are pathways.
#' @export
score_pathways <- function(data, mapping, method = c("mean_z", "median_z"), min_features = 3L, center = TRUE, scale = TRUE) {
  validate_omics_braid_data(data)
  method <- match.arg(method)
  req <- c("omic", "feature_id", "pathway")
  if (!all(req %in% names(mapping))) .stopf("mapping must contain columns: %s", paste(req, collapse = ", "))
  out <- list()
  coverage <- list()
  for (o in names(data$assays)) {
    mat <- data$assays[[o]]
    mp <- mapping[mapping$omic == o & mapping$feature_id %in% rownames(mat), req, drop = FALSE]
    if (!nrow(mp)) next
    splitf <- split(mp$feature_id, mp$pathway)
    splitf <- lapply(splitf, unique)
    splitf <- splitf[lengths(splitf) >= min_features]
    if (!length(splitf)) next
    z <- t(base::scale(t(mat), center = center, scale = scale))
    z[!is.finite(z)] <- NA_real_
    score <- matrix(NA_real_, nrow = length(splitf), ncol = ncol(mat), dimnames = list(names(splitf), colnames(mat)))
    covtab <- data.frame(omic = o, pathway = names(splitf), n_features = lengths(splitf), stringsAsFactors = FALSE)
    for (j in seq_along(splitf)) {
      xx <- z[splitf[[j]], , drop = FALSE]
      score[j, ] <- if (method == "mean_z") colMeans(xx, na.rm = TRUE) else apply(xx, 2L, stats::median, na.rm = TRUE)
      score[j, !is.finite(score[j, ])] <- NA_real_
    }
    out[[o]] <- score
    coverage[[o]] <- covtab
  }
  if (!length(out)) .stopf("No pathways met the mapping/min_features requirements.")
  ans <- omics_braid_data(out, data$metadata, sample_id = data$sample_id, annotation = do.call(rbind, coverage))
  attr(ans, "pathway_scoring") <- list(method = method, min_features = min_features, mapping_coverage = do.call(rbind, coverage))
  ans
}
