#' Empirically calibrate OmicsBraid omnibus and heterogeneity tests
#'
#' Provides resampling-based p-values for the two cross-omic quadratic tests.
#' The omnibus test can be calibrated by matched-subject label permutation or by
#' a centered matched-subject bootstrap. The heterogeneity test can be calibrated
#' by a raw-data null-shift matched bootstrap (recommended) or by an effect-level
#' centered bootstrap under the fitted common-effect null.
#'
#' The resampling is performed at the biological-subject level: all available
#' omic measurements belonging to a subject remain linked. This preserves the
#' cross-omic dependence that would be destroyed by shuffling individual assay
#' matrices independently.
#'
#' Empirical calibration is intended as a robust alternative when the chi-square
#' reference distributions used by `integrate_effects()` may be inaccurate, for
#' example under heavy-tailed sampling distributions. The asymptotic statistics
#' remain available and are not overwritten by this function.
#'
#' @param data An `omics_braid_data` object containing sample-level assays.
#' @param group Metadata column containing the two groups.
#' @param reference Reference-group label.
#' @param comparison Comparison-group label.
#' @param effects Optional layer-specific effect table. If supplied after
#'   `orient_omics()`, pass the same `orientation` so resampled effects receive
#'   the identical sign transformation.
#' @param entities Optional entities to calibrate. By default, entities observed
#'   in at least `min_omics` layers are used.
#' @param B Number of resampling replicates for each empirical null.
#' @param seed Random seed.
#' @param min_n Minimum observations per group within an omic.
#' @param min_omics Minimum omic layers per entity.
#' @param min_complete Minimum complete resampling draws required for a p-value.
#' @param omnibus_method Either `"permutation"` or `"centered_bootstrap"`.
#'   Permutation is appropriate for the global null in an exchangeable two-group
#'   design. Centered bootstrap is a nonparametric alternative.
#' @param heterogeneity_method `"null_shift_bootstrap"` (recommended robust
#'   calibration) or `"centered_bootstrap"`. Ordinary label permutation is not
#'   used because the heterogeneity null permits a common non-zero effect.
#' @param orientation Optional named +1/-1 vector applied to the resampled layer
#'   effects. This must match any scientific orientation already applied to
#'   `effects`.
#' @param p_adjust Multiple-testing method for empirical p-values.
#' @return A data frame containing asymptotic-independent empirical omnibus and
#'   heterogeneity p-values, empirical critical values, and resampling diagnostics.
#' @export
empirical_omics_tests <- function(data, group, reference, comparison,
                                  effects = NULL, entities = NULL,
                                  B = 499L, seed = 1L,
                                  min_n = 3L, min_omics = 2L,
                                  min_complete = 100L,
                                  omnibus_method = c("permutation", "centered_bootstrap"),
                                  heterogeneity_method = c("null_shift_bootstrap", "centered_bootstrap"),
                                  orientation = NULL,
                                  p_adjust = "BH") {
  validate_omics_braid_data(data)
  omnibus_method <- match.arg(omnibus_method)
  heterogeneity_method <- match.arg(heterogeneity_method)
  if (!group %in% names(data$metadata)) .stopf("Metadata column '%s' not found.", group)
  if (identical(as.character(reference), as.character(comparison))) .stopf("'reference' and 'comparison' must differ.")
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 99L) .stopf("'B' must be at least 99 for empirical calibration.")
  B <- as.integer(B)
  if (!is.numeric(min_complete) || length(min_complete) != 1L || min_complete < 50L) .stopf("'min_complete' must be at least 50.")
  if (min_complete > B) {
    min_complete <- max(50L, floor(B * .50))
    .warnf("'min_complete' exceeded B; using %d complete resamples as the minimum.", min_complete)
  }
  if (!is.null(orientation)) {
    if (!is.numeric(orientation) || is.null(names(orientation)) || any(!orientation %in% c(-1, 1))) {
      .stopf("'orientation' must be a named numeric vector containing only +1 and -1.")
    }
    if (anyDuplicated(names(orientation))) .stopf("'orientation' names must be unique.")
  }

  if (is.null(orientation) && !is.null(effects)) {
    orientation <- attr(effects, "orientation")
  }
  estimated_here <- is.null(effects)
  if (estimated_here) {
    effects <- estimate_effects(data, group, reference, comparison, entities = entities, min_n = min_n)
    if (!is.null(orientation)) effects <- orient_omics(effects, NULL, orientation)$effects
  } else {
    .validate_effect_table(effects)
  }

  counts <- table(effects$entity[is.finite(effects$effect) & is.finite(effects$se)])
  eligible <- names(counts[counts >= min_omics])
  if (is.null(entities)) entities <- eligible else entities <- intersect(as.character(entities), eligible)
  if (!length(entities)) .stopf("No entities have at least %d usable omics.", min_omics)

  # Generate only the resampling schemes that are needed. These draws are
  # independent of the matched bootstrap used by bootstrap_effect_covariance().
  perm_draws <- NULL
  boot_draws <- NULL
  nullq_draws <- NULL
  if (identical(omnibus_method, "permutation")) {
    perm_draws <- .resample_effect_draws(
      data, group, reference, comparison, entities = entities,
      B = B, seed = seed, min_n = min_n, mode = "permutation",
      orientation = orientation, return_z = TRUE
    )
  }
  # Ordinary stratified bootstrap is generated for the centered-bootstrap
  # omnibus comparator and, when requested, the centered-bootstrap Q null.
  boot_draws <- .resample_effect_draws(
    data, group, reference, comparison, entities = entities,
    B = B, seed = seed + 100003L, min_n = min_n,
    mode = "stratified_bootstrap", orientation = orientation,
    return_z = FALSE
  )
  if (identical(heterogeneity_method, "null_shift_bootstrap")) {
    null_data <- .impose_common_effect_null(
      data, group, reference, comparison, effects = effects,
      entities = entities, orientation = orientation, min_n = min_n
    )
    nullq_draws <- .resample_effect_draws(
      null_data, group, reference, comparison, entities = entities,
      B = B, seed = seed + 200003L, min_n = min_n,
      mode = "stratified_bootstrap", orientation = orientation,
      return_z = FALSE
    )
  }

  rows <- vector("list", length(entities)); kk <- 0L
  for (e in entities) {
    d <- effects[effects$entity == e & is.finite(effects$effect) & is.finite(effects$se) & effects$se > 0, , drop = FALSE]
    if (nrow(d) < min_omics) next
    os <- as.character(d$omic)
    y <- as.numeric(d$effect); names(y) <- os

    # Omnibus empirical null. When permutation is requested we also report a
    # centered-bootstrap calibration from the already generated heterogeneity
    # bootstrap draws, allowing direct robustness comparisons without a second
    # resampling pass.
    Xperm_z <- .collect_entity_resample_matrix(
      if (!is.null(perm_draws)) perm_draws$z else NULL, e, os
    )
    Xboot <- .collect_entity_resample_matrix(
      if (!is.null(boot_draws)) boot_draws$effect else NULL, e, os
    )
    zobs <- y / d$se
    op <- .permutation_omnibus_studentized(
      Xperm_z, zobs, min_complete = min_complete
    )
    ob <- .bootstrap_omnibus_studentized(
      Xboot, y, min_complete = min_complete
    )
    oo <- if (identical(omnibus_method, "permutation")) op else ob
    p_o <- oo$p
    crit_o <- oo$critical95
    stat_o <- oo$stat
    n_o <- oo$n
    cond_o <- oo$condition

    # Heterogeneity empirical null. Report both the effect-level centered
    # bootstrap and the stronger raw-data null-shift bootstrap when available.
    Xh_center <- .collect_entity_resample_matrix(
      if (!is.null(boot_draws)) boot_draws$effect else NULL, e, os
    )
    q_center <- .bootstrap_heterogeneity_deviation(
      Xh_center, y, min_complete = min_complete, recenter = TRUE
    )
    Xh_shift <- .collect_entity_resample_matrix(
      if (!is.null(nullq_draws)) nullq_draws$effect else NULL, e, os
    )
    q_shift <- .bootstrap_heterogeneity_deviation(
      Xh_shift, y, min_complete = min_complete, recenter = FALSE
    )
    qh <- if (identical(heterogeneity_method, "null_shift_bootstrap")) q_shift else q_center
    n_h <- qh$n
    p_h <- qh$p
    crit_h <- qh$critical95
    stat_h <- qh$stat
    theta_h <- qh$center
    cond_h <- qh$condition

    kk <- kk + 1L
    rows[[kk]] <- data.frame(
      entity = e,
      n_omics = length(os),
      omics = paste(os, collapse = ";"),
      empirical_omnibus_method = omnibus_method,
      empirical_heterogeneity_method = heterogeneity_method,
      T_omnibus_empirical = stat_o,
      p_omnibus_empirical = p_o,
      omnibus_empirical_critical_95 = crit_o,
      omnibus_resamples_used = n_o,
      omnibus_empirical_cov_condition = cond_o,
      p_omnibus_permutation = op$p,
      omnibus_permutation_critical_95 = op$critical95,
      omnibus_permutation_resamples_used = op$n,
      p_omnibus_centered_bootstrap = ob$p,
      omnibus_bootstrap_critical_95 = ob$critical95,
      omnibus_bootstrap_resamples_used = ob$n,
      common_effect_center = theta_h,
      T_heterogeneity_empirical = stat_h,
      p_heterogeneity_empirical = p_h,
      heterogeneity_empirical_critical_95 = crit_h,
      heterogeneity_resamples_used = n_h,
      heterogeneity_empirical_cov_condition = cond_h,
      p_heterogeneity_centered_bootstrap = q_center$p,
      p_heterogeneity_null_shift_bootstrap = q_shift$p,
      heterogeneity_centered_resamples_used = q_center$n,
      heterogeneity_null_shift_resamples_used = q_shift$n,
      stringsAsFactors = FALSE
    )
  }
  rows <- rows[seq_len(kk)]
  if (!length(rows)) .stopf("No entities had enough complete empirical resamples.")
  ans <- do.call(rbind, rows)
  ans$p_omnibus_empirical_adj <- stats::p.adjust(ans$p_omnibus_empirical, method = p_adjust)
  ans$p_heterogeneity_empirical_adj <- stats::p.adjust(ans$p_heterogeneity_empirical, method = p_adjust)
  rownames(ans) <- NULL
  attr(ans, "B") <- B
  attr(ans, "seed") <- seed
  attr(ans, "estimated_effects_internally") <- estimated_here
  ans
}

.quadratic_values <- function(X, W) {
  if (!nrow(X)) return(numeric())
  rowSums((X %*% W) * X)
}

.gls_heterogeneity_values <- function(X, W) {
  if (!nrow(X)) return(numeric())
  one <- rep(1, ncol(X))
  den <- as.numeric(t(one) %*% W %*% one)
  if (!is.finite(den) || den <= 0) return(rep(NA_real_, nrow(X)))
  theta <- as.numeric(X %*% W %*% one / den)
  R <- X - theta
  rowSums((R %*% W) * R)
}

.empirical_covariance_inverse <- function(X) {
  S <- try(stats::cov(X, use = "complete.obs"), silent = TRUE)
  if (inherits(S, "try-error") || any(!is.finite(S))) {
    return(list(W = NULL, condition = NA_real_))
  }
  S <- .make_psd(S)
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  cond <- max(ev) / max(min(ev), .Machine$double.eps)
  list(W = .safe_inverse(S), condition = cond)
}

.permutation_omnibus_studentized <- function(Z, zobs, min_complete) {
  if (is.null(Z) || !nrow(Z)) {
    return(list(p = NA_real_, critical95 = NA_real_, stat = NA_real_,
                n = 0L, condition = NA_real_))
  }
  Z <- Z[stats::complete.cases(Z), , drop = FALSE]
  if (nrow(Z) < min_complete || any(!is.finite(zobs))) {
    return(list(p = NA_real_, critical95 = NA_real_, stat = NA_real_,
                n = nrow(Z), condition = NA_real_))
  }
  # Under the joint null, group labels are exchangeable. The entire matched
  # subject is permuted, so cross-omic dependence is retained. We use a
  # Mahalanobis/Wald-type quadratic statistic with covariance estimated from
  # the joint permutation distribution. This is the empirical analogue of the
  # covariance-aware asymptotic omnibus statistic rather than a sum of
  # independent layer-wise z-squares.
  mu <- colMeans(Z)
  Zc <- sweep(Z, 2L, mu, "-")
  cw <- .empirical_covariance_inverse(Zc)
  if (is.null(cw$W)) {
    return(list(p = NA_real_, critical95 = NA_real_, stat = NA_real_,
                n = nrow(Z), condition = NA_real_))
  }
  zc <- zobs - mu
  stat <- as.numeric(t(zc) %*% cw$W %*% zc)
  tt <- .quadratic_values(Zc, cw$W)
  tt <- tt[is.finite(tt)]
  p <- if (length(tt) >= min_complete) (1 + sum(tt >= stat)) / (length(tt) + 1) else NA_real_
  critical <- if (length(tt) >= min_complete) {
    as.numeric(stats::quantile(tt, .95, names = FALSE, type = 8, na.rm = TRUE))
  } else NA_real_
  list(p = p, critical95 = critical, stat = stat, n = length(tt), condition = cw$condition)
}

.bootstrap_omnibus_studentized <- function(X, y, min_complete) {
  if (is.null(X) || !nrow(X)) {
    return(list(p = NA_real_, critical95 = NA_real_, stat = NA_real_,
                n = 0L, condition = NA_real_))
  }
  X <- X[stats::complete.cases(X), , drop = FALSE]
  if (nrow(X) < min_complete || any(!is.finite(y))) {
    return(list(p = NA_real_, critical95 = NA_real_, stat = NA_real_,
                n = nrow(X), condition = NA_real_))
  }
  # Center the ordinary matched bootstrap at its bootstrap mean to represent
  # sampling variation under the global zero-effect null. The covariance of
  # those centered effect draws supplies the empirical Wald metric.
  R <- sweep(X, 2L, colMeans(X), "-")
  cw <- .empirical_covariance_inverse(R)
  if (is.null(cw$W)) {
    return(list(p = NA_real_, critical95 = NA_real_, stat = NA_real_,
                n = nrow(X), condition = NA_real_))
  }
  stat <- as.numeric(t(y) %*% cw$W %*% y)
  tt <- .quadratic_values(R, cw$W)
  tt <- tt[is.finite(tt)]
  p <- if (length(tt) >= min_complete) (1 + sum(tt >= stat)) / (length(tt) + 1) else NA_real_
  critical <- if (length(tt) >= min_complete) {
    as.numeric(stats::quantile(tt, .95, names = FALSE, type = 8, na.rm = TRUE))
  } else NA_real_
  list(p = p, critical95 = critical, stat = stat, n = length(tt), condition = cw$condition)
}

.bootstrap_heterogeneity_deviation <- function(X, y, min_complete, recenter = TRUE) {
  if (is.null(X) || !nrow(X)) {
    return(list(p = NA_real_, critical95 = NA_real_, stat = NA_real_,
                n = 0L, condition = NA_real_, center = NA_real_))
  }
  X <- X[stats::complete.cases(X), , drop = FALSE]
  if (nrow(X) < min_complete || any(!is.finite(y))) {
    return(list(p = NA_real_, critical95 = NA_real_, stat = NA_real_,
                n = nrow(X), condition = NA_real_, center = NA_real_))
  }
  # The centered-bootstrap version removes the observed layer-specific means.
  # The null-shift version already comes from raw data transformed to have a
  # common standardized effect, so it is used without column recentering.
  Xnull <- if (isTRUE(recenter)) sweep(X, 2L, colMeans(X), "-") else X
  cw <- .empirical_covariance_inverse(Xnull)
  if (is.null(cw$W)) {
    return(list(p = NA_real_, critical95 = NA_real_, stat = NA_real_,
                n = nrow(X), condition = NA_real_, center = NA_real_))
  }
  one <- rep(1, length(y))
  den <- as.numeric(t(one) %*% cw$W %*% one)
  if (!is.finite(den) || den <= 0) {
    return(list(p = NA_real_, critical95 = NA_real_, stat = NA_real_,
                n = nrow(X), condition = cw$condition, center = NA_real_))
  }
  theta <- as.numeric(t(one) %*% cw$W %*% y / den)
  dy <- y - theta
  stat <- as.numeric(t(dy) %*% cw$W %*% dy)
  tt <- .gls_heterogeneity_values(Xnull, cw$W)
  tt <- tt[is.finite(tt)]
  p <- if (length(tt) >= min_complete) (1 + sum(tt >= stat)) / (length(tt) + 1) else NA_real_
  critical <- if (length(tt) >= min_complete) {
    as.numeric(stats::quantile(tt, .95, names = FALSE, type = 8, na.rm = TRUE))
  } else NA_real_
  list(p = p, critical95 = critical, stat = stat, n = length(tt), condition = cw$condition, center = theta)
}

.impose_common_effect_null <- function(data, group, reference, comparison,
                                       effects, entities, orientation = NULL,
                                       min_n = 3L) {
  validate_omics_braid_data(data)
  .validate_effect_table(effects)
  meta <- data$metadata
  ids <- as.character(meta[[data$sample_id]])
  grp_all <- as.character(meta[[group]])
  assays <- lapply(data$assays, function(x) x)

  theta <- vapply(entities, function(e) {
    d <- effects[effects$entity == e & is.finite(effects$effect) &
                   is.finite(effects$se) & effects$se > 0, , drop = FALSE]
    if (!nrow(d)) return(NA_real_)
    w <- 1 / (d$se^2)
    sum(w * d$effect) / sum(w)
  }, numeric(1))

  for (o in names(assays)) {
    mat <- assays[[o]]
    idx_meta <- match(colnames(mat), ids)
    grp <- grp_all[idx_meta]
    use_group <- grp %in% c(reference, comparison)
    i0 <- which(use_group & grp == reference)
    i1 <- which(use_group & grp == comparison)
    if (length(i0) < min_n || length(i1) < min_n) next
    feats <- intersect(entities, rownames(mat))
    mult <- if (!is.null(orientation) && o %in% names(orientation)) orientation[o] else 1
    for (e in feats) {
      th <- theta[e]
      if (!is.finite(th)) next
      x0 <- mat[e, i0]
      x1 <- mat[e, i1]
      x0 <- x0[is.finite(x0)]
      x1 <- x1[is.finite(x1)]
      n0 <- length(x0); n1 <- length(x1)
      if (n0 < min_n || n1 < min_n) next
      v0 <- stats::var(x0); v1 <- stats::var(x1)
      df <- n0 + n1 - 2
      sp2 <- ((n0 - 1) * v0 + (n1 - 1) * v1) / df
      if (!is.finite(sp2) || sp2 <= 0) next
      J <- .hedges_J(df)
      if (!is.finite(J) || J <= 0) next
      target_raw_g <- th / mult
      desired_diff <- (target_raw_g / J) * sqrt(sp2)
      observed_diff <- mean(x1) - mean(x0)
      delta <- desired_diff - observed_diff
      hit <- i1[is.finite(mat[e, i1])]
      mat[e, hit] <- mat[e, hit] + delta
    }
    assays[[o]] <- mat
  }
  omics_braid_data(assays = assays, metadata = meta, sample_id = data$sample_id)
}

.resample_effect_draws <- function(data, group, reference, comparison, entities,
                                   B, seed, min_n, mode,
                                   orientation = NULL, return_z = FALSE) {
  mode <- match.arg(mode, c("permutation", "stratified_bootstrap"))
  meta <- data$metadata
  ids <- as.character(meta[[data$sample_id]])
  grp <- as.character(meta[[group]])
  keep_meta <- grp %in% c(reference, comparison)
  ids <- ids[keep_meta]; grp <- grp[keep_meta]
  if (sum(grp == reference) < min_n || sum(grp == comparison) < min_n) {
    .stopf("Too few samples in one or both groups for empirical resampling.")
  }
  i_ref <- which(grp == reference)
  i_cmp <- which(grp == comparison)
  omics <- names(data$assays)
  out_effect <- stats::setNames(lapply(omics, function(o) {
    matrix(NA_real_, nrow = B, ncol = length(entities), dimnames = list(NULL, entities))
  }), omics)
  out_z <- if (isTRUE(return_z)) stats::setNames(lapply(omics, function(o) {
    matrix(NA_real_, nrow = B, ncol = length(entities), dimnames = list(NULL, entities))
  }), omics) else NULL
  set.seed(seed)
  for (b in seq_len(B)) {
    if (identical(mode, "permutation")) {
      bid <- ids
      bg <- sample(grp, length(grp), replace = FALSE)
    } else {
      ii <- c(sample(i_ref, length(i_ref), replace = TRUE),
              sample(i_cmp, length(i_cmp), replace = TRUE))
      bid <- ids[ii]
      bg <- grp[ii]
    }
    for (o in omics) {
      mat <- data$assays[[o]]
      present <- bid %in% colnames(mat)
      if (sum(present) < 2L * min_n) next
      sid <- bid[present]
      sg <- bg[present]
      if (sum(sg == reference) < min_n || sum(sg == comparison) < min_n) next
      feats <- intersect(entities, rownames(mat))
      if (!length(feats)) next
      idx <- match(sid, colnames(mat))
      bm <- mat[feats, idx, drop = FALSE]
      gz <- .hedges_g_se_matrix(bm, sg, reference, comparison, min_n = min_n)
      gv <- gz$g
      zv <- gz$z
      if (!is.null(orientation) && o %in% names(orientation)) {
        gv <- gv * orientation[o]
        zv <- zv * orientation[o]
      }
      out_effect[[o]][b, feats] <- gv
      if (isTRUE(return_z)) out_z[[o]][b, feats] <- zv
    }
  }
  list(effect = out_effect, z = out_z)
}

.collect_entity_resample_matrix <- function(draws, entity, omics) {
  if (is.null(draws)) return(matrix(NA_real_, nrow = 0L, ncol = length(omics), dimnames = list(NULL, omics)))
  B <- max(vapply(draws, nrow, integer(1)))
  X <- matrix(NA_real_, nrow = B, ncol = length(omics), dimnames = list(NULL, omics))
  for (j in seq_along(omics)) {
    o <- omics[j]
    M <- draws[[o]]
    if (!is.null(M) && entity %in% colnames(M)) X[seq_len(nrow(M)), j] <- M[, entity]
  }
  X
}

.apply_empirical_test_results <- function(integrated, empirical,
                                          use_empirical_as_primary = FALSE,
                                          p_adjust = "BH") {
  out <- integrated
  ii <- match(out$entity, empirical$entity)
  out$p_omnibus_asymptotic <- out$p_omnibus
  out$p_heterogeneity_asymptotic <- out$p_heterogeneity
  out$p_omnibus_empirical <- empirical$p_omnibus_empirical[ii]
  out$p_heterogeneity_empirical <- empirical$p_heterogeneity_empirical[ii]
  out$empirical_omnibus_method <- empirical$empirical_omnibus_method[ii]
  out$empirical_heterogeneity_method <- empirical$empirical_heterogeneity_method[ii]
  out$omnibus_resamples_used <- empirical$omnibus_resamples_used[ii]
  out$heterogeneity_resamples_used <- empirical$heterogeneity_resamples_used[ii]
  if (isTRUE(use_empirical_as_primary)) {
    use_o <- is.finite(out$p_omnibus_empirical)
    use_h <- is.finite(out$p_heterogeneity_empirical)
    out$p_omnibus[use_o] <- out$p_omnibus_empirical[use_o]
    out$p_heterogeneity[use_h] <- out$p_heterogeneity_empirical[use_h]
    out$p_omnibus_adj <- stats::p.adjust(out$p_omnibus, method = p_adjust)
    out$p_heterogeneity_adj <- stats::p.adjust(out$p_heterogeneity, method = p_adjust)
    out$heterogeneity_test_flag <- is.finite(out$p_heterogeneity_adj) & out$p_heterogeneity_adj <= 0.05
    out$has_omnibus_evidence_local <- is.finite(out$p_omnibus) & out$p_omnibus <= 0.05
    out$has_omnibus_evidence <- is.finite(out$p_omnibus_adj) & out$p_omnibus_adj <= 0.05
    out$evidence_qualified_direction_agreement <- ifelse(
      out$has_omnibus_evidence, out$direction_agreement, NA_real_
    )
    out$test_p_method <- ifelse(use_o | use_h, "empirical_when_available", "asymptotic")
  } else {
    out$test_p_method <- "asymptotic_primary_empirical_reported"
  }
  out
}
