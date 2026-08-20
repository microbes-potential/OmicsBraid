#' Integrate effects across omic layers and quantify heterogeneity
#'
#' Uses generalized least squares (GLS) to estimate a common cross-omic effect.
#' A generalized Cochran Q statistic tests whether the layer-specific effects are
#' compatible with a common effect after accounting for their sampling covariance.
#' A separate multivariate Wald-type omnibus statistic tests the joint null that
#' all layer effects are zero; unlike the consensus effect, this test does not
#' cancel equally strong effects occurring in opposite directions.
#' The reported I2-like statistic is descriptive and should not be interpreted as
#' literal between-study heterogeneity because omics layers are not studies.
#'
#' @param effects Data frame containing `entity`, `omic`, `effect`, and `se`.
#' @param covariance Optional output of `bootstrap_effect_covariance()` or a named list of covariance matrices.
#' @param min_omics Minimum omics per entity.
#' @param p_adjust Multiple-testing method for integrated and heterogeneity p-values.
#' @param conf_level Confidence level for the analytic GLS consensus interval.
#' @return Data frame of integrated effects and heterogeneity diagnostics.
#' @export
integrate_effects <- function(effects, covariance = NULL, min_omics = 2L, p_adjust = "BH", conf_level = 0.95) {
  .validate_effect_table(effects)
  if (!is.numeric(conf_level) || length(conf_level) != 1L || !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) .stopf("'conf_level' must lie strictly between 0 and 1.")
  crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  covlist <- NULL
  covmode <- "independence"
  if (!is.null(covariance)) {
    if (inherits(covariance, "omics_braid_covariance")) {
      covlist <- covariance$covariance
      covmode <- "matched_bootstrap"
    } else if (is.list(covariance)) {
      covlist <- covariance
      covmode <- "supplied"
    } else .stopf("Unsupported covariance object.")
  }
  ents <- unique(effects$entity)
  out <- vector("list", length(ents)); k <- 0L
  for (e in ents) {
    d <- effects[effects$entity == e & is.finite(effects$effect) & is.finite(effects$se) & effects$se > 0, , drop = FALSE]
    if (nrow(d) < min_omics) next
    y <- d$effect; names(y) <- d$omic
    V <- diag(d$se^2, nrow = nrow(d), ncol = nrow(d))
    dimnames(V) <- list(d$omic, d$omic)
    entity_covmode <- "independence"
    if (!is.null(covlist) && !is.null(covlist[[e]])) {
      V0 <- covlist[[e]]
      os <- intersect(d$omic, rownames(V0))
      if (length(os)) {
        V[os, os] <- V0[os, os, drop = FALSE]
        entity_covmode <- if (length(os) == nrow(d)) covmode else paste0("partial_", covmode)
      }
    }
    V <- .make_psd(V)
    W <- .safe_inverse(V)
    one <- rep(1, length(y))
    den <- as.numeric(t(one) %*% W %*% one)
    if (!is.finite(den) || den <= 0) next
    gls_weights <- as.numeric(W %*% one / den)
    theta <- sum(gls_weights * y)
    se <- sqrt(1 / den)
    z <- theta / se
    p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
    W_omnibus <- as.numeric(t(y) %*% W %*% y)
    df_omnibus <- length(y)
    p_omnibus <- stats::pchisq(W_omnibus, df = df_omnibus, lower.tail = FALSE)
    Q <- as.numeric(t(y - theta) %*% W %*% (y - theta))
    df <- length(y) - 1L
    phet <- stats::pchisq(Q, df = df, lower.tail = FALSE)
    I2 <- if (Q > 0) max(0, (Q - df) / Q) * 100 else 0
    wmarg <- 1 / (d$se^2)
    if (abs(theta) < .Machine$double.eps) {
      diragree <- max(sum(wmarg[y >= 0]), sum(wmarg[y <= 0])) / sum(wmarg)
    } else {
      diragree <- sum(wmarg[sign(y) == sign(theta)]) / sum(wmarg)
    }
    ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
    cond_v <- max(ev) / max(min(ev), .Machine$double.eps)
    k <- k + 1L
    out[[k]] <- data.frame(entity = e, n_omics = length(y), omics = paste(d$omic, collapse = ";"),
                           integrated_effect = theta, integrated_se = se,
                           conf_low = theta - crit * se,
                           conf_high = theta + crit * se,
                           z_value = z, p_value = p,
                           W_omnibus = W_omnibus, df_omnibus = df_omnibus,
                           p_omnibus = p_omnibus,
                           Q_omics = Q, df_heterogeneity = df,
                           p_heterogeneity = phet, I2_omics = I2,
                           direction_agreement = diragree,
                           min_gls_weight = min(gls_weights), max_gls_weight = max(gls_weights),
                           has_negative_gls_weight = any(gls_weights < -1e-8),
                           covariance_condition_number = cond_v,
                           covariance_mode = entity_covmode,
                           stringsAsFactors = FALSE)
  }
  out <- out[seq_len(k)]
  if (!length(out)) .stopf("No entities had at least %d usable omics.", min_omics)
  ans <- do.call(rbind, out)
  ans$p_adj <- stats::p.adjust(ans$p_value, method = p_adjust)
  ans$p_omnibus_adj <- stats::p.adjust(ans$p_omnibus, method = p_adjust)
  ans$p_heterogeneity_adj <- stats::p.adjust(ans$p_heterogeneity, method = p_adjust)
  ans$heterogeneity_test_flag <- is.finite(ans$p_heterogeneity_adj) & ans$p_heterogeneity_adj < 0.05
  ans$has_omnibus_evidence_local <- is.finite(ans$p_omnibus) & ans$p_omnibus < 0.05
  ans$has_omnibus_evidence <- is.finite(ans$p_omnibus_adj) & ans$p_omnibus_adj < 0.05
  # Raw directional alignment can be high by chance under the null. This
  # evidence-qualified version is intentionally NA when the joint-null test is
  # not rejected, preventing null/noisy entities from appearing perfectly
  # concordant in overview graphics.
  ans$evidence_qualified_direction_agreement <- ifelse(
    ans$has_omnibus_evidence, ans$direction_agreement, NA_real_
  )
  rownames(ans) <- NULL
  ans
}
