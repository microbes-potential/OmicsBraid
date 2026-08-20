#' Harmonize assay-specific feature identifiers to common entities
#'
#' Maps feature IDs within each omic to a shared entity identifier for
#' entity-level cross-omic analysis (for example Ensembl RNA identifiers and
#' UniProt proteins to a common gene symbol). Many-to-one mappings are rejected
#' by default because collapsing isoforms/probes changes the scientific estimand.
#'
#' @param data An `omics_braid_data` object.
#' @param mapping Data frame with columns `omic`, `feature_id`, and `entity`.
#' @param collapse How to handle multiple assay features mapping to one entity:
#'   `"error"` (default), `"mean"`, or `"median"`.
#' @param min_mapped Minimum mapped entities required per retained omic.
#' @return An `omics_braid_data` object with harmonized entity row names and a
#'   mapping report stored in `attr(x, "entity_harmonization")`.
#' @export
harmonize_entities <- function(data, mapping,
                               collapse = c("error", "mean", "median"),
                               min_mapped = 1L) {
  validate_omics_braid_data(data)
  collapse <- match.arg(collapse)
  req <- c("omic", "feature_id", "entity")
  if (!is.data.frame(mapping) || !all(req %in% names(mapping))) {
    .stopf("'mapping' must be a data.frame containing: %s", paste(req, collapse = ", "))
  }
  mapping$omic <- as.character(mapping$omic)
  mapping$feature_id <- as.character(mapping$feature_id)
  mapping$entity <- as.character(mapping$entity)
  if (anyNA(mapping[, req]) || any(!nzchar(mapping$omic)) || any(!nzchar(mapping$feature_id)) || any(!nzchar(mapping$entity))) {
    .stopf("Mapping omic, feature_id, and entity values must be non-missing and non-empty.")
  }
  if (anyDuplicated(mapping[c("omic", "feature_id")])) {
    .stopf("Each omic + feature_id may map to only one common entity in one harmonization call.")
  }

  out <- list(); report <- list()
  for (o in names(data$assays)) {
    mat <- data$assays[[o]]
    mp <- mapping[mapping$omic == o & mapping$feature_id %in% rownames(mat), req, drop = FALSE]
    if (!nrow(mp)) next
    mp <- mp[match(intersect(rownames(mat), mp$feature_id), mp$feature_id), , drop = FALSE]
    mat <- mat[mp$feature_id, , drop = FALSE]
    target <- mp$entity
    dup_target <- duplicated(target) | duplicated(target, fromLast = TRUE)
    if (any(dup_target) && collapse == "error") {
      ex <- unique(target[dup_target])
      .stopf("Omic '%s' has multiple features mapping to common entities (e.g. %s). Choose collapse='mean'/'median' only if scientifically justified.",
             o, paste(utils::head(ex, 5L), collapse = ", "))
    }
    if (any(dup_target)) {
      lev <- unique(target)
      cm <- matrix(NA_real_, nrow = length(lev), ncol = ncol(mat), dimnames = list(lev, colnames(mat)))
      for (j in seq_along(lev)) {
        xx <- mat[target == lev[j], , drop = FALSE]
        if (nrow(xx) == 1L) cm[j, ] <- xx[1L, ]
        else if (collapse == "mean") cm[j, ] <- colMeans(xx, na.rm = TRUE)
        else cm[j, ] <- apply(xx, 2L, stats::median, na.rm = TRUE)
      }
      cm[!is.finite(cm)] <- NA_real_
      mat2 <- cm
    } else {
      rownames(mat) <- target
      mat2 <- mat
    }
    if (nrow(mat2) < min_mapped) next
    out[[o]] <- mat2
    report[[o]] <- data.frame(
      omic = o,
      input_features = nrow(data$assays[[o]]),
      mapped_features = nrow(mp),
      mapped_entities = nrow(mat2),
      collapsed = sum(dup_target),
      stringsAsFactors = FALSE
    )
  }
  if (length(out) < 1L) .stopf("No omic retained at least %d mapped entities.", min_mapped)
  ans <- omics_braid_data(out, data$metadata, sample_id = data$sample_id,
                          annotation = mapping)
  attr(ans, "entity_harmonization") <- do.call(rbind, report)
  ans
}
