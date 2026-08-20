#' Test ordered cross-omic effect trajectories with generalized least squares
#'
#' Fits a covariance-aware linear trajectory to standardized effects across an
#' explicitly ordered set of omic layers. Effects are aligned to the dominant
#' observed direction before fitting, so a positive slope represents increasing
#' absolute effect magnitude (amplification) and a negative slope represents
#' decreasing magnitude (attenuation). Three practical hypotheses are evaluated:
#' a meaningfully positive slope, a meaningfully negative slope, and practical
#' equivalence of the slope to a flat trajectory within
#' `[-trajectory_margin, +trajectory_margin]`.
#'
#' The trend test is intended for ordered layers when attenuation/amplification
#' is scientifically meaningful. It does not establish causality or temporal
#' direction.
#'
#' @param effects Data frame containing `entity`, `omic`, `effect`, and `se`.
#' @param covariance Optional output of `bootstrap_effect_covariance()` or a named
#'   list of entity-specific covariance matrices. If omitted, layer estimates are
#'   treated as independent for this calculation.
#' @param omic_order Ordered character vector describing the layer trajectory.
#' @param trajectory_margin Smallest meaningful change in standardized effect per
#'   one-layer transition. A scalar greater than zero.
#' @param alpha Local significance level used to define the trend state.
#' @param min_omics Minimum observed layers required.
#' @param p_adjust Multiple-testing method used for confirmatory adjusted trend
#'   p-values across entities. Local states remain the default for braid geometry.
#' @return One row per entity containing the aligned GLS slope, uncertainty,
#'   practical trend tests, and local/adjusted trajectory states.
#' @export
test_braid_trend <- function(effects, covariance = NULL, omic_order,
                             trajectory_margin = 0.15, alpha = 0.05,
                             min_omics = 2L, p_adjust = "BH") {
  .validate_effect_table(effects)
  if (!is.character(omic_order) || length(omic_order) < 2L || anyDuplicated(omic_order)) {
    .stopf("'omic_order' must contain at least two unique omic names.")
  }
  if (!is.numeric(trajectory_margin) || length(trajectory_margin) != 1L ||
      !is.finite(trajectory_margin) || trajectory_margin <= 0) {
    .stopf("'trajectory_margin' must be one finite number > 0.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    .stopf("'alpha' must lie strictly between 0 and 1.")
  }

  covlist <- NULL
  covmode <- "independence"
  if (!is.null(covariance)) {
    if (inherits(covariance, "omics_braid_covariance")) {
      covlist <- covariance$covariance
      covmode <- "matched_bootstrap"
    } else if (is.list(covariance)) {
      covlist <- covariance
      covmode <- "supplied"
    } else {
      .stopf("Unsupported covariance object.")
    }
  }

  ents <- unique(effects$entity)
  out <- list(); kk <- 0L
  for (e in ents) {
    d <- effects[
      effects$entity == e & effects$omic %in% omic_order &
        is.finite(effects$effect) & is.finite(effects$se) & effects$se > 0,
      , drop = FALSE
    ]
    oo <- intersect(omic_order, d$omic)
    d <- d[match(oo, d$omic), , drop = FALSE]
    if (nrow(d) < min_omics) next

    y <- d$effect
    names(y) <- d$omic
    V <- diag(d$se^2, nrow = nrow(d), ncol = nrow(d))
    dimnames(V) <- list(d$omic, d$omic)
    entity_covmode <- "independence"
    if (!is.null(covlist) && !is.null(covlist[[e]])) {
      V0 <- covlist[[e]]
      if (!is.null(rownames(V0)) && !is.null(colnames(V0))) {
        os <- intersect(d$omic, intersect(rownames(V0), colnames(V0)))
        if (length(os)) {
          V[os, os] <- V0[os, os, drop = FALSE]
          entity_covmode <- if (length(os) == nrow(d)) covmode else paste0("partial_", covmode)
        }
      }
    }
    V <- .make_psd(V)
    W <- .safe_inverse(V)

    # Align all effects to the dominant precision-weighted direction. This makes
    # the slope describe magnitude propagation for both positive and negative
    # concordant trajectories. Classification only treats this slope as
    # confirmatory when the layer states themselves support a common direction.
    precision <- 1 / (d$se^2)
    direction <- sign(sum(precision * y, na.rm = TRUE))
    if (!is.finite(direction) || direction == 0) direction <- sign(mean(y, na.rm = TRUE))
    if (!is.finite(direction) || direction == 0) direction <- 1
    ya <- direction * y

    layer_index <- seq_len(nrow(d)) - 1
    X <- cbind(intercept = 1, layer = layer_index)
    XtWX <- t(X) %*% W %*% X
    Bcov <- .safe_inverse(XtWX)
    beta <- as.numeric(Bcov %*% t(X) %*% W %*% ya)
    slope <- beta[2]
    slope_se <- sqrt(max(Bcov[2, 2], 0))
    if (!is.finite(slope_se) || slope_se <= 0) next

    z0 <- slope / slope_se
    p0 <- 2 * stats::pnorm(abs(z0), lower.tail = FALSE)

    # One-sided tests against a scientifically meaningful trajectory margin.
    z_amp <- (slope - trajectory_margin) / slope_se
    p_amp <- stats::pnorm(z_amp, lower.tail = FALSE) # H0: slope <= +margin
    z_att <- (slope + trajectory_margin) / slope_se
    p_att <- stats::pnorm(z_att, lower.tail = TRUE)  # H0: slope >= -margin

    # TOST for a practically flat trajectory.
    z_flat_lower <- (slope + trajectory_margin) / slope_se
    z_flat_upper <- (slope - trajectory_margin) / slope_se
    p_flat_lower <- stats::pnorm(z_flat_lower, lower.tail = FALSE)
    p_flat_upper <- stats::pnorm(z_flat_upper, lower.tail = TRUE)
    p_flat <- max(p_flat_lower, p_flat_upper)

    kk <- kk + 1L
    out[[kk]] <- data.frame(
      entity = e,
      n_layers_trend = nrow(d),
      trend_omics = paste(d$omic, collapse = ";"),
      aligned_direction = as.integer(direction),
      same_sign_estimates = all(y > 0) || all(y < 0),
      trend_intercept = beta[1],
      trend_slope = slope,
      trend_slope_se = slope_se,
      trend_conf_low = slope - stats::qnorm(.975) * slope_se,
      trend_conf_high = slope + stats::qnorm(.975) * slope_se,
      p_trend_zero = p0,
      p_meaningful_amplification = p_amp,
      p_meaningful_attenuation = p_att,
      p_flat_tost = p_flat,
      trajectory_margin = trajectory_margin,
      covariance_mode_trend = entity_covmode,
      stringsAsFactors = FALSE
    )
  }

  if (!length(out)) {
    return(data.frame(
      entity = character(), n_layers_trend = integer(), trend_omics = character(),
      aligned_direction = integer(), same_sign_estimates = logical(),
      trend_intercept = numeric(), trend_slope = numeric(), trend_slope_se = numeric(),
      trend_conf_low = numeric(), trend_conf_high = numeric(), p_trend_zero = numeric(),
      p_meaningful_amplification = numeric(), p_meaningful_attenuation = numeric(),
      p_flat_tost = numeric(), trajectory_margin = numeric(),
      covariance_mode_trend = character(), trend_state_local = character(),
      p_meaningful_amplification_adj = numeric(), p_meaningful_attenuation_adj = numeric(),
      p_flat_tost_adj = numeric(), trend_state_adjusted = character(),
      stringsAsFactors = FALSE
    ))
  }

  ans <- do.call(rbind, out)
  trend_state <- function(pamp, patt, pflat) {
    st <- rep("uncertain", length(pamp))
    st[is.finite(pflat) & pflat < alpha] <- "flat"
    st[st == "uncertain" & is.finite(pamp) & pamp < alpha] <- "amplification"
    st[st == "uncertain" & is.finite(patt) & patt < alpha] <- "attenuation"
    st
  }
  ans$trend_state_local <- trend_state(
    ans$p_meaningful_amplification,
    ans$p_meaningful_attenuation,
    ans$p_flat_tost
  )
  # A magnitude trajectory is not interpretable when estimated layer effects
  # already cross zero in opposite directions; inversion/discordance inference
  # should be used instead.
  ans$trend_state_local[!ans$same_sign_estimates] <- "not_applicable"
  ans$p_meaningful_amplification_adj <- stats::p.adjust(ans$p_meaningful_amplification, method = p_adjust)
  ans$p_meaningful_attenuation_adj <- stats::p.adjust(ans$p_meaningful_attenuation, method = p_adjust)
  ans$p_flat_tost_adj <- stats::p.adjust(ans$p_flat_tost, method = p_adjust)
  ans$trend_state_adjusted <- trend_state(
    ans$p_meaningful_amplification_adj,
    ans$p_meaningful_attenuation_adj,
    ans$p_flat_tost_adj
  )
  ans$trend_state_adjusted[!ans$same_sign_estimates] <- "not_applicable"
  rownames(ans) <- NULL
  ans
}
