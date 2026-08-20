#' Construct an OmicsBraid data object
#'
#' @param assays Named list of numeric matrices. Features are rows and samples are columns.
#' @param metadata Data frame containing one row per biological sample.
#' @param sample_id Column in `metadata` holding sample identifiers.
#' @param annotation Optional feature annotation data frame.
#' @return An object of class `omics_braid_data`.
#' @export
omics_braid_data <- function(assays, metadata, sample_id = "sample_id", annotation = NULL) {
  if (!is.list(assays) || !length(assays)) .stopf("'assays' must be a non-empty named list.")
  if (is.null(names(assays)) || any(!nzchar(names(assays))) || anyDuplicated(names(assays))) .stopf("'assays' must have unique non-empty names.")
  assay_names <- names(assays)
  assays <- lapply(seq_along(assays), function(i) .as_numeric_matrix(assays[[i]], assay_names[i]))
  names(assays) <- assay_names
  if (!is.data.frame(metadata)) metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  if (!sample_id %in% names(metadata)) .stopf("metadata does not contain sample ID column '%s'.", sample_id)
  ids <- as.character(metadata[[sample_id]])
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) .stopf("metadata sample IDs must be non-missing and unique.")
  metadata[[sample_id]] <- ids
  obj <- structure(list(assays = assays, metadata = metadata, sample_id = sample_id, annotation = annotation), class = "omics_braid_data")
  validate_omics_braid_data(obj)
  obj
}

#' Coerce supported objects to OmicsBraid data
#' @param x Object to coerce.
#' @param ... Additional arguments.
#' @export
as_omics_braid_data <- function(x, ...) UseMethod("as_omics_braid_data")

#' @export
as_omics_braid_data.list <- function(x, metadata, sample_id = "sample_id", annotation = NULL, ...) {
  omics_braid_data(x, metadata = metadata, sample_id = sample_id, annotation = annotation)
}

#' @export
as_omics_braid_data.MultiAssayExperiment <- function(x, sample_id = "sample_id", annotation = NULL,
                                                       replicate_action = c("error", "mean"), ...) {
  if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) .stopf("Install Bioconductor package 'MultiAssayExperiment' for this coercion method.")
  replicate_action <- match.arg(replicate_action)
  aa <- MultiAssayExperiment::assays(x)
  assays <- lapply(aa, as.matrix)
  sm <- as.data.frame(MultiAssayExperiment::sampleMap(x), stringsAsFactors = FALSE)
  cd <- as.data.frame(MultiAssayExperiment::colData(x), stringsAsFactors = FALSE)
  if (!sample_id %in% names(cd)) cd[[sample_id]] <- rownames(cd)

  for (o in names(assays)) {
    mat <- assays[[o]]
    mp <- sm[as.character(sm$assay) == o, , drop = FALSE]
    primary <- as.character(mp$primary[match(colnames(mat), as.character(mp$colname))])
    if (anyNA(primary)) {
      bad <- colnames(mat)[is.na(primary)]
      .warnf("Dropping %d columns in assay '%s' that are absent from the MultiAssayExperiment sampleMap.", length(bad), o)
      mat <- mat[, !is.na(primary), drop = FALSE]
      primary <- primary[!is.na(primary)]
    }
    if (anyDuplicated(primary)) {
      if (replicate_action == "error") {
        .stopf("Assay '%s' has multiple assay columns mapping to the same primary subject. Re-run with replicate_action='mean' only if averaging those replicates is scientifically appropriate.", o)
      }
      lev <- unique(primary)
      collapsed <- matrix(NA_real_, nrow(mat), length(lev), dimnames = list(rownames(mat), lev))
      for (j in seq_along(lev)) {
        xx <- mat[, primary == lev[j], drop = FALSE]
        collapsed[, j] <- if (ncol(xx) == 1L) xx[,1L] else rowMeans(xx, na.rm = TRUE)
      }
      collapsed[!is.finite(collapsed)] <- NA_real_
      mat <- collapsed
    } else {
      colnames(mat) <- primary
    }
    assays[[o]] <- mat
  }
  omics_braid_data(assays, metadata = cd, sample_id = sample_id, annotation = annotation)
}

#' Validate an OmicsBraid data object
#' @param x An `omics_braid_data` object.
#' @return Invisibly returns `TRUE` or throws an informative error.
#' @export
validate_omics_braid_data <- function(x) {
  if (!inherits(x, "omics_braid_data")) .stopf("Object is not of class 'omics_braid_data'.")
  ids <- as.character(x$metadata[[x$sample_id]])
  overlap <- vapply(x$assays, function(a) sum(colnames(a) %in% ids), numeric(1))
  if (any(overlap == 0L)) .stopf("At least one assay has no sample IDs matching metadata: %s", paste(names(overlap)[overlap == 0L], collapse = ", "))
  extra <- lapply(x$assays, function(a) setdiff(colnames(a), ids))
  if (any(lengths(extra) > 0L)) {
    bad <- names(extra)[lengths(extra) > 0L]
    .warnf("Some assay columns are absent from metadata and will be ignored: %s", paste(bad, collapse = ", "))
  }
  invisible(TRUE)
}

#' @export
print.omics_braid_data <- function(x, ...) {
  cat("<omics_braid_data>\n")
  cat(" Samples in metadata:", nrow(x$metadata), "\n")
  cat(" Assays:", paste(names(x$assays), collapse = ", "), "\n")
  for (nm in names(x$assays)) cat(sprintf("  - %s: %d features x %d samples\n", nm, nrow(x$assays[[nm]]), ncol(x$assays[[nm]])))
  invisible(x)
}
