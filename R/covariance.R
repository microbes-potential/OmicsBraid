#' Estimate cross-omic effect covariance by matched-subject bootstrap
#'
#' Biological subjects are resampled as whole units, preserving their matched
#' measurements across omics. Bootstrap correlations are estimated for each
#' entity and combined with analytic marginal standard errors from
#' `estimate_effects()`: V = diag(SE) R_boot diag(SE). This stabilizes marginal
#' uncertainty while retaining empirically estimated cross-layer dependence.
#'
#' @param data An `omics_braid_data` object.
#' @param group Metadata group column.
#' @param reference Reference group.
#' @param comparison Comparison group.
#' @param effects Optional output of `estimate_effects()`.
#' @param entities Optional entities to bootstrap. By default, entities present in at least two omics.
#' @param B Number of subject-level bootstrap replicates.
#' @param seed Random seed.
#' @param min_n Minimum observations per group within an omic.
#' @param shrinkage Correlation shrinkage toward the identity in [0,1].
#' @param min_complete Minimum usable bootstrap pairs to estimate a correlation.
#' @param stratified Logical; resample subjects separately within reference and comparison groups.
#'   This preserves the observed group sizes and is recommended for fixed two-group designs.
#' @return Object of class `omics_braid_covariance` containing a covariance matrix per entity.
#' @export
bootstrap_effect_covariance <- function(data, group, reference, comparison, effects = NULL,
                                        entities = NULL, B = 500L, seed = 1L, min_n = 3L,
                                        shrinkage = 0.05, min_complete = 50L,
                                        stratified = TRUE) {
  validate_omics_braid_data(data)
  if (!group %in% names(data$metadata)) .stopf("Metadata column '%s' not found.", group)
  if (is.null(effects)) effects <- estimate_effects(data, group, reference, comparison, min_n = min_n)
  if (B < 50L) .warnf("B < 50 gives unstable covariance estimates; >= 300 is recommended for final analysis.")
  if (shrinkage < 0 || shrinkage > 1) .stopf("'shrinkage' must lie in [0,1].")
  counts <- table(effects$entity[is.finite(effects$effect)])
  eligible <- names(counts[counts >= 2L])
  if (is.null(entities)) entities <- eligible else entities <- intersect(entities, eligible)
  if (!length(entities)) .stopf("No entities are observed in at least two omics.")
  if (length(entities) > 2000L) .warnf("Bootstrapping %d entities may be computationally expensive. Pathway-level or prefiltered analysis is recommended.", length(entities))
  meta <- data$metadata
  ids <- as.character(meta[[data$sample_id]])
  grp_all <- as.character(meta[[group]])
  use_meta <- grp_all %in% c(reference, comparison)
  ids <- ids[use_meta]; grp_all <- grp_all[use_meta]
  n <- length(ids)
  if (n < 2L * min_n) .stopf("Too few metadata samples for bootstrap.")
  set.seed(seed)
  omics <- names(data$assays)
  boot <- stats::setNames(lapply(omics, function(o) matrix(NA_real_, nrow = B, ncol = length(entities), dimnames = list(NULL, entities))), omics)
  i_ref <- which(grp_all == reference)
  i_cmp <- which(grp_all == comparison)
  for (b in seq_len(B)) {
    if (isTRUE(stratified)) {
      ii <- c(sample(i_ref, length(i_ref), replace = TRUE),
              sample(i_cmp, length(i_cmp), replace = TRUE))
    } else {
      ii <- sample.int(n, n, replace = TRUE)
    }
    bid <- ids[ii]; bg <- grp_all[ii]
    for (o in omics) {
      mat <- data$assays[[o]]
      keep <- bid %in% colnames(mat)
      if (sum(keep) < 2L * min_n) next
      sid <- bid[keep]; sg <- bg[keep]
      if (sum(sg == reference) < min_n || sum(sg == comparison) < min_n) next
      feats <- intersect(entities, rownames(mat))
      if (!length(feats)) next
      bm <- mat[feats, match(sid, colnames(mat)), drop = FALSE]
      gv <- .hedges_g_matrix(bm, sg, reference, comparison, min_n = min_n)
      boot[[o]][b, feats] <- gv
    }
  }
  covs <- list(); cors <- list(); n_pairs <- list(); diagnostics <- list()
  for (e in entities) {
    ef <- effects[effects$entity == e & is.finite(effects$effect) & is.finite(effects$se), , drop = FALSE]
    os <- intersect(ef$omic, omics)
    if (length(os) < 2L) next
    X <- sapply(os, function(o) boot[[o]][, e])
    if (is.null(dim(X))) X <- matrix(X, ncol = length(os))
    colnames(X) <- os
    R <- diag(length(os)); dimnames(R) <- list(os, os)
    NP <- matrix(0L, length(os), length(os), dimnames = list(os, os))
    for (i in seq_along(os)) for (j in seq_along(os)) {
      ok <- is.finite(X[,i]) & is.finite(X[,j])
      NP[i,j] <- sum(ok)
      if (i != j && sum(ok) >= min_complete) {
        rr <- suppressWarnings(stats::cor(X[ok,i], X[ok,j]))
        if (is.finite(rr)) R[i,j] <- rr else R[i,j] <- 0
      }
    }
    R <- (1 - shrinkage) * R + shrinkage * diag(length(os))
    dimnames(R) <- list(os, os)
    R_before_psd <- R
    R <- .make_psd(R)
    psd_delta <- sqrt(sum((R - R_before_psd)^2))
    se <- ef$se[match(os, ef$omic)]
    V <- outer(se, se) * R
    dimnames(V) <- list(os, os)
    covs[[e]] <- .make_psd(V)
    cors[[e]] <- R
    n_pairs[[e]] <- NP
    offdiag <- R[row(R) != col(R)]
    pair_counts <- NP[upper.tri(NP)]
    diagnostics[[e]] <- data.frame(
      entity = e, n_omics = length(os),
      min_boot_pairs = min(pair_counts, na.rm = TRUE),
      n_pairs_below_min_complete = sum(pair_counts < min_complete, na.rm = TRUE),
      max_abs_correlation = if (length(offdiag)) max(abs(offdiag), na.rm = TRUE) else 0,
      psd_adjustment_frobenius = psd_delta,
      stringsAsFactors = FALSE
    )
  }
  diagtab <- if (length(diagnostics)) do.call(rbind, diagnostics) else data.frame()
  structure(list(covariance = covs, correlation = cors, n_pairs = n_pairs,
                 diagnostics = diagtab, B = B, seed = seed, shrinkage = shrinkage,
                 min_complete = min_complete, min_n = min_n, stratified = stratified, entities = names(covs),
                 boot_effects = boot), class = "omics_braid_covariance")
}
