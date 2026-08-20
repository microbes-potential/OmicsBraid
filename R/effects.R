#' Estimate layer-specific standardized effects
#'
#' Estimates Hedges' g for a two-group contrast in each feature/pathway and omic.
#' Inputs should already be quality-controlled and normalized appropriately for
#' their assay technology. Hedges' g is scale-free but does not repair poor raw
#' preprocessing, severe censoring, or inappropriate transformations.
#'
#' @param data An `omics_braid_data` object.
#' @param group Metadata column containing the two groups.
#' @param reference Reference group label.
#' @param comparison Comparison group label. Positive effects mean comparison > reference.
#' @param entities Optional character vector limiting features/pathways.
#' @param min_n Minimum non-missing observations per group and omic.
#' @param conf_level Confidence level.
#' @param p_adjust Multiple-testing method applied separately within each omic.
#' @return Data frame of effect estimates and uncertainty.
#' @export
estimate_effects <- function(data, group, reference, comparison, entities = NULL, min_n = 3L,
                             conf_level = 0.95, p_adjust = "BH") {
  validate_omics_braid_data(data)
  if (!group %in% names(data$metadata)) .stopf("Metadata column '%s' not found.", group)
  if (identical(as.character(reference), as.character(comparison))) .stopf("'reference' and 'comparison' must be different groups.")
  if (!is.numeric(conf_level) || length(conf_level) != 1L || !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) .stopf("'conf_level' must lie strictly between 0 and 1.")
  if (!is.numeric(min_n) || length(min_n) != 1L || min_n < 2L) .stopf("'min_n' must be at least 2.")
  ids <- as.character(data$metadata[[data$sample_id]])
  grp_all <- as.character(data$metadata[[group]])
  if (!reference %in% grp_all || !comparison %in% grp_all) .stopf("Both reference and comparison groups must occur in metadata.")
  alpha <- 1 - conf_level
  out <- list()
  for (o in names(data$assays)) {
    mat <- data$assays[[o]]
    keep_samples <- colnames(mat) %in% ids
    mat <- mat[, keep_samples, drop = FALSE]
    meta_idx <- match(colnames(mat), ids)
    grp <- grp_all[meta_idx]
    keep_group <- grp %in% c(reference, comparison)
    mat <- mat[, keep_group, drop = FALSE]
    grp <- grp[keep_group]
    if (!is.null(entities)) mat <- mat[intersect(rownames(mat), entities), , drop = FALSE]
    if (!nrow(mat)) next
    i0 <- which(grp == reference); i1 <- which(grp == comparison)
    x0 <- mat[, i0, drop = FALSE]; x1 <- mat[, i1, drop = FALSE]
    n0 <- rowSums(!is.na(x0)); n1 <- rowSums(!is.na(x1))
    m0 <- .row_mean_na(x0); m1 <- .row_mean_na(x1)
    v0 <- .row_var_na(x0); v1 <- .row_var_na(x1)
    df <- n0 + n1 - 2
    sp2 <- ((n0 - 1) * v0 + (n1 - 1) * v1) / df
    d <- (m1 - m0) / sqrt(sp2)
    J <- .hedges_J(df)
    g <- J * d
    # Common large-sample approximation to Var(g); the bootstrap covariance
    # function uses these analytic marginal SEs and estimates dependence separately.
    varg <- (n0 + n1) / (n0 * n1) + (g * g) / (2 * pmax(n0 + n1, 1))
    se <- sqrt(varg)
    z <- g / se
    p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
    crit <- stats::qnorm(1 - alpha / 2)
    low <- g - crit * se; high <- g + crit * se
    bad <- n0 < min_n | n1 < min_n | !is.finite(g) | !is.finite(se) | sp2 <= 0
    g[bad] <- se[bad] <- varg[bad] <- p[bad] <- low[bad] <- high[bad] <- NA_real_
    tab <- data.frame(entity = rownames(mat), omic = o, effect = g, se = se, var = varg,
                      conf_low = low, conf_high = high, p_value = p, n_reference = n0,
                      n_comparison = n1, df = df, reference = reference,
                      comparison = comparison, stringsAsFactors = FALSE)
    out[[o]] <- tab
  }
  if (!length(out)) .stopf("No effects could be estimated.")
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans$p_adj <- .adjust_within_omic(ans$p_value, ans$omic, method = p_adjust)
  attr(ans, "effect_measure") <- "Hedges_g"
  attr(ans, "conf_level") <- conf_level
  ans
}
