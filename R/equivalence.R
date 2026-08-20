#' Test practical equivalence to a negligible effect region
#'
#' Performs two one-sided tests (TOST) for whether each standardized effect lies
#' within `[-margin, +margin]`. It reports both local (entity-specific) and
#' multiplicity-adjusted inferential states. Local states are recommended for
#' describing the geometry of an individual braid because they do not change
#' merely when unrelated entities are added to the analysis; adjusted states are
#' retained for confirmatory screening across many entities.
#'
#' @param effects Effect table from `estimate_effects()` or compatible summary statistics.
#' @param margin Smallest effect size of interest on the standardized-effect scale; scalar or named by omic.
#' @param alpha Significance level.
#' @param p_adjust Multiple-testing adjustment applied separately within each omic.
#' @param state_basis Which inferential state is copied to the legacy `state`
#'   column: `"local"` (default) or `"adjusted"`.
#' @return Effect table with TOST p-values, local/adjusted difference p-values,
#'   and both local and multiplicity-adjusted practical states.
#' @export
test_equivalence <- function(effects, margin = 0.3, alpha = 0.05,
                             p_adjust = "BH",
                             state_basis = c("local", "adjusted")) {
  state_basis <- match.arg(state_basis)
  .validate_effect_table(effects)
  m <- .resolve_margin(margin, effects$omic)
  est <- effects$effect; se <- effects$se
  lo <- -m; hi <- m
  zl <- (est - lo) / se
  zu <- (est - hi) / se
  # Large-sample normal approximation, consistent with the standardized-effect
  # inference used elsewhere in the current research prototype.
  p_lower <- stats::pnorm(zl, lower.tail = FALSE)
  p_upper <- stats::pnorm(zu, lower.tail = TRUE)
  p_tost <- pmax(p_lower, p_upper)
  p_tost[!is.finite(est) | !is.finite(se) | se <= 0] <- NA_real_
  p_tost_adj <- .adjust_within_omic(p_tost, effects$omic, method = p_adjust)

  p_diff <- if ("p_value" %in% names(effects)) effects$p_value else {
    2 * stats::pnorm(abs(est / se), lower.tail = FALSE)
  }
  p_diff_adj <- if ("p_adj" %in% names(effects)) effects$p_adj else {
    .adjust_within_omic(p_diff, effects$omic, method = p_adjust)
  }

  make_state <- function(eq_p, diff_p) {
    equivalent <- is.finite(eq_p) & eq_p < alpha
    different <- is.finite(diff_p) & diff_p < alpha
    state <- rep("uncertain", length(est))
    state[equivalent] <- "equivalent"
    state[!equivalent & different & est > 0] <- "positive"
    state[!equivalent & different & est < 0] <- "negative"
    list(state = state, equivalent = equivalent)
  }

  local <- make_state(p_tost, p_diff)
  adjusted <- make_state(p_tost_adj, p_diff_adj)

  out <- effects
  out$equiv_margin <- m
  out$p_tost_lower <- p_lower
  out$p_tost_upper <- p_upper
  out$p_tost <- p_tost
  out$p_tost_adj <- p_tost_adj
  out$p_difference <- p_diff
  out$p_difference_adj <- p_diff_adj
  out$equivalent_local <- local$equivalent
  out$equivalent_adjusted <- adjusted$equivalent
  out$state_local <- local$state
  out$state_adjusted <- adjusted$state
  # Backward-compatible columns. The default is deliberately local so the
  # braid label for an entity is not altered by the number of unrelated tests.
  out$equivalent <- if (state_basis == "local") local$equivalent else adjusted$equivalent
  out$state <- if (state_basis == "local") local$state else adjusted$state
  out$state_basis <- state_basis
  out
}
