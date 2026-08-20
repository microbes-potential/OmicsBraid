#' Classify cross-omic braid patterns using inferential evidence
#'
#' Braid labels are deliberately conservative. Opposite statistically supported
#' layer directions confirm inversion. Buffering and emergence require practical
#' equivalence in the appropriate downstream/upstream layers. Concordance,
#' attenuation, and amplification require all observed layers to support one
#' direction and use a covariance-aware GLS trajectory test. If the joint-null
#' omnibus test is not rejected and practical equivalence is not established,
#' the result is labelled `no_detectable_effect` rather than incorrectly claiming
#' equivalence.
#'
#' A separate `suggestive_pattern` is derived from effect geometry only. It can
#' be useful when the inferential pattern is unresolved because the data are too
#' imprecise, but it is not a confirmatory conclusion.
#'
#' @param equivalence Output of `test_equivalence()`.
#' @param omic_order Ordered character vector describing the biological/display order of omics.
#' @param covariance Optional covariance object used when a trend table must be
#'   computed internally.
#' @param integrated Optional output of `integrate_effects()`. Supplying it allows
#'   `no_detectable_effect` to be distinguished from generic uncertainty.
#' @param trend Optional output of `test_braid_trend()`. If absent it is computed.
#' @param trajectory_margin Smallest meaningful effect change per one-layer
#'   transition for attenuation/amplification.
#' @param alpha Local inferential significance level.
#' @param min_slope Deprecated alias for `trajectory_margin` retained for early
#'   OmicsBraid prototypes.
#' @return One row per entity containing confirmatory and suggestive labels plus
#'   trajectory diagnostics.
#' @export
classify_braids <- function(equivalence, omic_order, covariance = NULL,
                            integrated = NULL, trend = NULL,
                            trajectory_margin = 0.15, alpha = 0.05,
                            min_slope = NULL) {
  req <- c("entity", "omic", "effect", "se", "state", "equiv_margin")
  if (!all(req %in% names(equivalence))) .stopf("equivalence table is missing required columns.")
  .validate_effect_table(equivalence)
  if (anyDuplicated(omic_order)) .stopf("'omic_order' must not contain duplicates.")
  if (!is.null(min_slope)) {
    trajectory_margin <- min_slope
    .warnf("'min_slope' is deprecated; use 'trajectory_margin'.")
  }
  if (!is.numeric(trajectory_margin) || length(trajectory_margin) != 1L ||
      !is.finite(trajectory_margin) || trajectory_margin <= 0) {
    .stopf("'trajectory_margin' must be one finite number > 0.")
  }
  if (is.null(trend)) {
    trend <- test_braid_trend(
      equivalence, covariance = covariance, omic_order = omic_order,
      trajectory_margin = trajectory_margin, alpha = alpha
    )
  }

  ents <- unique(equivalence$entity)
  out <- vector("list", length(ents))
  for (ii in seq_along(ents)) {
    e <- ents[ii]
    d <- equivalence[equivalence$entity == e & equivalence$omic %in% omic_order, , drop = FALSE]
    oo <- intersect(omic_order, d$omic)
    d <- d[match(oo, d$omic), , drop = FALSE]
    states <- as.character(d$state)
    vals <- d$effect
    margins <- d$equiv_margin

    suggestive <- if (nrow(d) >= 2L) {
      .pattern_from_values(vals, margins = margins, trajectory_margin = trajectory_margin)
    } else {
      "insufficient"
    }

    tr <- if (!is.null(trend) && nrow(trend)) trend[trend$entity == e, , drop = FALSE] else data.frame()
    omnibus_p <- NA_real_
    omnibus_evidence <- NA
    if (!is.null(integrated) && is.data.frame(integrated) && "entity" %in% names(integrated)) {
      ir <- integrated[integrated$entity == e, , drop = FALSE]
      if (nrow(ir) && "p_omnibus" %in% names(ir)) {
        omnibus_p <- ir$p_omnibus[1]
        omnibus_evidence <- is.finite(omnibus_p) && omnibus_p <= alpha
      }
    }

    pat <- "uncertain"
    status <- "unresolved"
    basis <- "insufficient inferential evidence"

    if (nrow(d) < 2L) {
      pat <- "insufficient"
      status <- "insufficient"
      basis <- "fewer than two observed layers"
    } else {
      sig_idx <- which(states %in% c("positive", "negative"))
      sigstates <- states[sig_idx]
      has_pos <- any(sigstates == "positive")
      has_neg <- any(sigstates == "negative")
      all_eq <- length(states) > 0L && all(states == "equivalent")

      if (all_eq) {
        pat <- "null_equivalent"
        status <- "confirmed"
        basis <- "all observed layers passed practical-equivalence testing"
      } else if (has_pos && has_neg) {
        pat <- "inversion"
        status <- "confirmed"
        basis <- "statistically supported effects occur in opposite directions"
      } else if (length(sig_idx) == 0L && identical(omnibus_evidence, FALSE)) {
        pat <- "no_detectable_effect"
        status <- "no_evidence"
        basis <- "joint-null omnibus test not rejected; equivalence not established"
      } else if (length(sig_idx) > 0L && length(unique(sigstates)) == 1L) {
        # Confirmed buffering: at least one directional layer followed only by
        # practically equivalent downstream layers.
        last_sig <- max(sig_idx)
        first_sig <- min(sig_idx)
        buffer_ok <- last_sig < nrow(d) && all(states[seq.int(last_sig + 1L, nrow(d))] == "equivalent")
        emerge_ok <- first_sig > 1L && all(states[seq_len(first_sig - 1L)] == "equivalent")

        if (buffer_ok && !emerge_ok) {
          pat <- "buffering"
          status <- "confirmed"
          basis <- "directional upstream effect followed by confirmed practical equivalence"
        } else if (emerge_ok && !buffer_ok) {
          pat <- "emergence"
          status <- "confirmed"
          basis <- "confirmed practical equivalence precedes a directional downstream effect"
        } else if (!buffer_ok && !emerge_ok && all(states %in% sigstates[1])) {
          # Every observed layer supports the same direction. Trend inference now
          # determines whether magnitude is flat, attenuating, or amplifying.
          tstate <- if (nrow(tr) && "trend_state_local" %in% names(tr)) tr$trend_state_local[1] else "uncertain"
          if (identical(tstate, "attenuation")) {
            pat <- "attenuation"
            status <- "confirmed"
            basis <- "same-direction layers plus meaningful negative GLS trajectory"
          } else if (identical(tstate, "amplification")) {
            pat <- "amplification"
            status <- "confirmed"
            basis <- "same-direction layers plus meaningful positive GLS trajectory"
          } else if (identical(tstate, "flat")) {
            pat <- if (sigstates[1] == "positive") "concordant_increase" else "concordant_decrease"
            status <- "confirmed"
            basis <- "same-direction layers plus GLS trajectory practically equivalent to flat"
          } else {
            # Concordance here means statistically supported direction across all
            # observed layers, not equality of magnitudes. If a meaningful GLS
            # trajectory cannot be established, retain the broader concordant
            # label and explicitly mark its trajectory subtype as unresolved.
            pat <- if (sigstates[1] == "positive") "concordant_increase" else "concordant_decrease"
            status <- "direction_confirmed"
            basis <- "all observed layers support one direction; magnitude trajectory is unresolved"
          }
        } else {
          pat <- "uncertain"
          status <- "unresolved"
          basis <- "partial or transient layer pattern does not satisfy a confirmatory braid rule"
        }
      } else if (length(sig_idx) == 0L && isTRUE(omnibus_evidence)) {
        pat <- "uncertain"
        status <- "unresolved"
        basis <- "joint multi-omic evidence exists but no individual layer state is resolved"
      }
    }

    interpretation <- if (status == "confirmed") {
      pat
    } else if (status == "direction_confirmed") {
      paste0(pat, " (trajectory unresolved; geometry: ", suggestive, ")")
    } else if (status == "no_evidence") {
      if (!suggestive %in% c("uncertain", "insufficient")) {
        paste0("no_detectable_effect (geometry: ", suggestive, ")")
      } else {
        "no_detectable_effect"
      }
    } else if (!suggestive %in% c("uncertain", "insufficient")) {
      paste0("suggestive_", suggestive)
    } else {
      pat
    }

    row <- data.frame(
      entity = e,
      pattern = pat,
      pattern_status = status,
      suggestive_pattern = suggestive,
      interpretation_label = interpretation,
      confirmatory_basis = basis,
      n_layers_observed = nrow(d),
      n_positive = sum(states == "positive"),
      n_negative = sum(states == "negative"),
      n_equivalent = sum(states == "equivalent"),
      n_uncertain = sum(states == "uncertain"),
      omnibus_p_local = omnibus_p,
      omnibus_evidence_local = omnibus_evidence,
      stringsAsFactors = FALSE
    )

    if (nrow(tr)) {
      keep <- c(
        "trend_slope", "trend_slope_se", "trend_conf_low", "trend_conf_high",
        "p_trend_zero", "p_meaningful_amplification", "p_meaningful_attenuation",
        "p_flat_tost", "trajectory_margin", "trend_state_local",
        "p_meaningful_amplification_adj", "p_meaningful_attenuation_adj",
        "p_flat_tost_adj", "trend_state_adjusted", "covariance_mode_trend"
      )
      keep <- intersect(keep, names(tr))
      for (nm in keep) row[[nm]] <- tr[[nm]][1]
    }
    out[[ii]] <- row
  }
  do.call(rbind, out)
}

#' Quantify uncertainty in effect-braid geometry
#'
#' Draws from a multivariate normal approximation to the estimated cross-omic
#' effects and reports how often each practical geometric pattern occurs. These
#' frequencies propagate estimation uncertainty; they are not Bayesian posterior
#' probabilities and they do not replace the confirmatory inferential braid label.
#'
#' @param effects Effect table.
#' @param covariance Optional bootstrap covariance object/list. If absent, omics are treated as independent.
#' @param omic_order Ordered omics.
#' @param margin Equivalence/negligible-effect margin; scalar or named by omic.
#' @param n_draws Number of Monte Carlo draws per entity.
#' @param trajectory_margin Practical trajectory threshold per layer transition.
#' @param seed Random seed.
#' @param min_slope Deprecated alias for `trajectory_margin`.
#' @return Long data frame of geometric pattern frequencies per entity.
#' @export
braid_pattern_probabilities <- function(effects, covariance = NULL, omic_order,
                                        margin = 0.3, n_draws = 2000L,
                                        trajectory_margin = 0.15, seed = 1L,
                                        min_slope = NULL) {
  .validate_effect_table(effects)
  if (anyDuplicated(omic_order)) .stopf("'omic_order' must not contain duplicates.")
  if (!is.null(min_slope)) {
    trajectory_margin <- min_slope
    .warnf("'min_slope' is deprecated; use 'trajectory_margin'.")
  }
  covlist <- NULL
  if (inherits(covariance, "omics_braid_covariance")) covlist <- covariance$covariance else if (is.list(covariance)) covlist <- covariance
  set.seed(seed)
  ents <- unique(effects$entity)
  out <- list(); z <- 0L
  for (e in ents) {
    d <- effects[effects$entity == e & effects$omic %in% omic_order & is.finite(effects$effect) & is.finite(effects$se), , drop = FALSE]
    oo <- intersect(omic_order, d$omic)
    d <- d[match(oo, d$omic), , drop = FALSE]
    if (nrow(d) < 2L) next
    mu <- d$effect; names(mu) <- d$omic
    if (!is.null(covlist) && !is.null(covlist[[e]]) &&
        !is.null(rownames(covlist[[e]])) && all(d$omic %in% rownames(covlist[[e]]))) {
      V <- covlist[[e]][d$omic, d$omic, drop = FALSE]
    } else {
      V <- diag(d$se^2, nrow = nrow(d), ncol = nrow(d))
      dimnames(V) <- list(d$omic, d$omic)
    }
    mm <- .resolve_margin(margin, d$omic)
    X <- .rmvn_psd(n_draws, mu, V)
    pats <- apply(X, 1L, .pattern_from_values,
                  margins = mm, trajectory_margin = trajectory_margin)
    tt <- sort(table(pats), decreasing = TRUE)
    z <- z + 1L
    out[[z]] <- data.frame(
      entity = e, pattern = names(tt), count = as.integer(tt),
      probability = as.numeric(tt) / n_draws, stringsAsFactors = FALSE
    )
  }
  if (!length(out)) return(data.frame(entity = character(), pattern = character(), count = integer(), probability = numeric()))
  do.call(rbind, out)
}

# Summarize geometric uncertainty draws relative to inferential/suggestive labels.
.summarize_pattern_uncertainty <- function(classification, probabilities) {
  if (is.null(probabilities) || !nrow(probabilities)) return(classification)
  pieces <- lapply(seq_len(nrow(classification)), function(i) {
    row <- classification[i, , drop = FALSE]
    pp <- probabilities[probabilities$entity == row$entity, , drop = FALSE]
    if (!nrow(pp)) {
      row$uncertainty_mode_pattern <- NA_character_
      row$uncertainty_mode_support <- NA_real_
      row$deterministic_pattern_support <- NA_real_
      row$suggestive_pattern_support <- NA_real_
      row$pattern_entropy <- NA_real_
      row$mode_agrees_with_deterministic <- NA
      row$classification_stability <- if (row$pattern_status %in% c("unresolved", "insufficient")) "unresolved" else "not_estimated"
      return(row)
    }
    pp <- pp[order(pp$probability, decreasing = TRUE), , drop = FALSE]
    mode_pat <- pp$pattern[1]
    mode_prob <- pp$probability[1]
    hit_det <- pp$probability[pp$pattern == row$pattern]
    det_support <- if (length(hit_det)) hit_det[1] else NA_real_
    hit_sug <- pp$probability[pp$pattern == row$suggestive_pattern]
    sug_support <- if (length(hit_sug)) hit_sug[1] else 0
    pr <- pp$probability[is.finite(pp$probability) & pp$probability > 0]
    ent <- if (length(pr) <= 1L) 0 else -sum(pr * log(pr)) / log(length(pr))
    agree <- if (row$pattern_status %in% c("confirmed", "direction_confirmed")) identical(as.character(mode_pat), as.character(row$pattern)) else NA
    stability <- if (row$pattern_status %in% c("unresolved", "insufficient")) {
      "unresolved"
    } else if (row$pattern_status == "no_evidence") {
      if (sug_support >= 0.80) "geometry_high" else if (sug_support >= 0.60) "geometry_moderate" else "geometry_low"
    } else if (isTRUE(agree) && is.finite(det_support) && det_support >= 0.80) {
      if (row$pattern_status == "direction_confirmed") "direction_high" else "high"
    } else if (isTRUE(agree) && is.finite(det_support) && det_support >= 0.60) {
      if (row$pattern_status == "direction_confirmed") "direction_moderate" else "moderate"
    } else {
      if (row$pattern_status == "direction_confirmed") "direction_low" else "low"
    }
    row$uncertainty_mode_pattern <- mode_pat
    row$uncertainty_mode_support <- mode_prob
    row$deterministic_pattern_support <- det_support
    row$suggestive_pattern_support <- sug_support
    row$pattern_entropy <- ent
    row$mode_agrees_with_deterministic <- agree
    row$classification_stability <- stability
    row
  })
  do.call(rbind, pieces)
}
