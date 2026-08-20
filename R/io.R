#' Read OmicsBraid input files
#'
#' @param metadata_file CSV/TSV file containing sample metadata.
#' @param assay_files Named character vector or named list mapping omic names to CSV/TSV files.
#' @param sample_id Metadata sample identifier column.
#' @param sep Separator. If `NULL`, inferred from file extension (`.tsv`/`.txt` -> tab; otherwise comma).
#' @param annotation_file Optional annotation CSV/TSV file.
#' @return An `omics_braid_data` object.
#' @export
read_omics_braid <- function(metadata_file, assay_files, sample_id = "sample_id", sep = NULL, annotation_file = NULL) {
  read_one <- function(path) {
    s <- sep
    if (is.null(s)) s <- if (grepl("\\.(tsv|txt)$", path, ignore.case = TRUE)) "\t" else ","
    utils::read.table(path, header = TRUE, sep = s, check.names = FALSE, stringsAsFactors = FALSE, quote = "\"", comment.char = "")
  }
  meta <- read_one(metadata_file)
  if (is.null(names(assay_files)) || any(!nzchar(names(assay_files)))) .stopf("'assay_files' must be named by omic.")
  assays <- lapply(assay_files, function(path) {
    tab <- read_one(path)
    if (ncol(tab) < 2L) .stopf("Assay file '%s' must contain a feature ID column plus sample columns.", path)
    rn <- as.character(tab[[1L]])
    mat <- as.matrix(tab[-1L])
    storage.mode(mat) <- "double"
    rownames(mat) <- rn
    mat
  })
  ann <- if (!is.null(annotation_file)) read_one(annotation_file) else NULL
  omics_braid_data(assays, meta, sample_id = sample_id, annotation = ann)
}
