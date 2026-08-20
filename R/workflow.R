#' Run the complete OmicsBraid workflow on sample-level data
#'
#' @param data An `omics_braid_data` object.
#' @param group Metadata group column.
#' @param reference Reference group.
#' @param comparison Comparison group.
#' @param omic_order Biological/display order of omics. Defaults to assay order.
#' @param entities Optional entities to analyze.
#' @param pathway_mapping Optional mapping passed to `score_pathways()`; if supplied, analysis is performed on pathway scores.
#' @param pathway_method Pathway scoring method.
#' @param min_pathway_features Minimum pathway features per omic.
#' @param bootstrap_B Bootstrap replicates for cross-omic covariance. Set to 0 to assume independence.
#' @param bootstrap_shrinkage Correlation shrinkage.
#' @param empirical_tests Logical; if `TRUE`, calculate additional resampling-
#'   calibrated omnibus and heterogeneity p-values using `empirical_omics_tests()`.
#' @param empirical_B Number of empirical resampling replicates.
#' @param empirical_omnibus_method `"permutation"` or `"centered_bootstrap"`.
#' @param empirical_heterogeneity_method `"null_shift_bootstrap"` (recommended)
#'   or `"centered_bootstrap"`.
#' @param empirical_use_as_primary Logical; if `TRUE`, empirical p-values replace
#'   asymptotic p-values when available for downstream omnibus evidence. The
#'   original asymptotic p-values are retained in separate columns.
#' @param ci_method Layer-effect confidence interval method: `"analytic"`,
#'   `"percentile"`, `"basic"`, or `"bca"`. Bootstrap intervals change interval
#'   reporting only; analytic SEs and p-values remain available.
#' @param integrated_ci_method Consensus-effect confidence interval method: `"analytic"`, `"percentile"`, or `"basic"`.
#' @param ci_conf_level Confidence level for analytic/bootstrap intervals.
#' @param ci_min_boot Minimum finite bootstrap draws required before a bootstrap interval replaces the analytic interval.
#' @param orientation Optional named +1/-1 vector to harmonize omic effect directions.
#' @param equivalence_margin Smallest effect size of interest for practical equivalence.
#' @param alpha Significance level for local difference/equivalence and braid inference.
#' @param state_basis Use local or multiplicity-adjusted inferential states for deterministic braid labels.
#' @param equivalence_p_adjust Multiple-testing method retained for adjusted equivalence/difference evidence.
#' @param trajectory_margin Smallest meaningful standardized-effect change per one-layer transition.
#' @param trend_p_adjust Multiple-testing method retained for adjusted trajectory evidence.
#' @param min_slope Deprecated alias for `trajectory_margin`.
#' @param pattern_draws Monte Carlo draws for geometric uncertainty propagation.
#' @param seed Random seed.
#' @return An object of class `omics_braid_result`.
#' @export
run_omics_braid <- function(data, group, reference, comparison,
                            omic_order = names(data$assays), entities = NULL,
                            pathway_mapping = NULL, pathway_method = "mean_z",
                            min_pathway_features = 3L, bootstrap_B = 500L,
                            bootstrap_shrinkage = 0.05,
                            empirical_tests = FALSE, empirical_B = 499L,
                            empirical_omnibus_method = c("permutation", "centered_bootstrap"),
                            empirical_heterogeneity_method = c("null_shift_bootstrap", "centered_bootstrap"),
                            empirical_use_as_primary = FALSE,
                            ci_method = c("analytic", "percentile", "basic", "bca"),
                            integrated_ci_method = c("analytic", "percentile", "basic"),
                            ci_conf_level = 0.95, ci_min_boot = 100L,
                            orientation = NULL,
                            equivalence_margin = 0.3, alpha = 0.05,
                            state_basis = c("local", "adjusted"), equivalence_p_adjust = "BH",
                            trajectory_margin = 0.15, trend_p_adjust = "BH",
                            min_slope = NULL, pattern_draws = 2000L, seed = 1L) {
  state_basis <- match.arg(state_basis)
  ci_method <- match.arg(ci_method)
  integrated_ci_method <- match.arg(integrated_ci_method)
  empirical_omnibus_method <- match.arg(empirical_omnibus_method)
  empirical_heterogeneity_method <- match.arg(empirical_heterogeneity_method)
  if (!is.numeric(ci_conf_level) || length(ci_conf_level) != 1L || !is.finite(ci_conf_level) || ci_conf_level <= 0 || ci_conf_level >= 1) {
    .stopf("'ci_conf_level' must lie strictly between 0 and 1.")
  }
  if (ci_method != "analytic" && bootstrap_B <= 0L) .stopf("Bootstrap layer CIs require 'bootstrap_B' > 0.")
  if (integrated_ci_method != "analytic" && bootstrap_B <= 0L) .stopf("Bootstrap consensus CIs require 'bootstrap_B' > 0.")
  if (!is.null(min_slope)) {
    trajectory_margin <- min_slope
    .warnf("'min_slope' is deprecated; use 'trajectory_margin'.")
  }
  dat <- data
  analysis_level <- "feature"
  if (!is.null(pathway_mapping)) {
    dat <- score_pathways(data, pathway_mapping, method = pathway_method, min_features = min_pathway_features)
    analysis_level <- "pathway"
  }
  eff <- estimate_effects(dat, group, reference, comparison, entities = entities, conf_level = ci_conf_level)
  covobj <- NULL
  if (bootstrap_B > 0L) {
    covobj <- bootstrap_effect_covariance(dat, group, reference, comparison, effects = eff,
                                          entities = entities, B = bootstrap_B, seed = seed,
                                          shrinkage = bootstrap_shrinkage)
  }
  if (!is.null(orientation)) {
    oriented <- orient_omics(eff, covobj, orientation)
    eff <- oriented$effects; covobj <- oriented$covariance
  }

  effect_intervals <- NULL
  if (ci_method != "analytic") {
    effect_intervals <- bootstrap_effect_intervals(
      effects = eff, bootstrap = covobj, method = ci_method,
      conf_level = ci_conf_level, min_boot = ci_min_boot,
      data = if (ci_method == "bca") dat else NULL,
      group = if (ci_method == "bca") group else NULL,
      reference = if (ci_method == "bca") reference else NULL,
      comparison = if (ci_method == "bca") comparison else NULL
    )
    eff <- .apply_layer_bootstrap_intervals(eff, effect_intervals, ci_method)
  } else {
    eff$ci_method <- "analytic"
  }

  integ <- integrate_effects(eff, covariance = covobj, conf_level = ci_conf_level)

  empirical <- NULL
  if (isTRUE(empirical_tests)) {
    empirical <- empirical_omics_tests(
      data = dat, group = group, reference = reference, comparison = comparison,
      effects = eff, entities = entities, B = empirical_B, seed = seed + 1009L,
      min_complete = max(50L, floor(empirical_B * .50)),
      omnibus_method = empirical_omnibus_method,
      heterogeneity_method = empirical_heterogeneity_method,
      orientation = orientation
    )
    integ <- .apply_empirical_test_results(
      integ, empirical,
      use_empirical_as_primary = empirical_use_as_primary,
      p_adjust = "BH"
    )
  }

  consensus_intervals <- NULL
  if (integrated_ci_method != "analytic") {
    consensus_intervals <- bootstrap_consensus_intervals(
      effects = eff, bootstrap = covobj, integrated = integ,
      method = integrated_ci_method, conf_level = ci_conf_level,
      min_boot = ci_min_boot
    )
    integ <- .apply_consensus_bootstrap_intervals(integ, consensus_intervals, integrated_ci_method)
  } else {
    integ$ci_method <- "analytic"
  }

  eq <- test_equivalence(eff, margin = equivalence_margin, alpha = alpha,
                         p_adjust = equivalence_p_adjust, state_basis = state_basis)
  tr <- test_braid_trend(
    eff, covariance = covobj, omic_order = omic_order,
    trajectory_margin = trajectory_margin, alpha = alpha,
    p_adjust = trend_p_adjust
  )
  cls <- classify_braids(
    eq, omic_order = omic_order, covariance = covobj,
    integrated = integ, trend = tr,
    trajectory_margin = trajectory_margin, alpha = alpha
  )
  probs <- braid_pattern_probabilities(
    eff, covariance = covobj, omic_order = omic_order,
    margin = equivalence_margin, n_draws = pattern_draws,
    trajectory_margin = trajectory_margin, seed = seed + 1L
  )
  cls2 <- .summarize_pattern_uncertainty(cls, probs)
  out <- structure(list(data = dat, effects = eff, covariance = covobj,
                        effect_intervals = effect_intervals,
                        integrated = integ, empirical_tests = empirical,
                        consensus_intervals = consensus_intervals,
                        equivalence = eq, trend = tr,
                        classification = cls2, pattern_probabilities = probs,
                        settings = list(group = group, reference = reference,
                                        comparison = comparison, omic_order = omic_order,
                                        equivalence_margin = equivalence_margin,
                                        alpha = alpha, state_basis = state_basis,
                                        equivalence_p_adjust = equivalence_p_adjust,
                                        trajectory_margin = trajectory_margin,
                                        trend_p_adjust = trend_p_adjust,
                                        bootstrap_B = bootstrap_B,
                                        bootstrap_shrinkage = bootstrap_shrinkage,
                                        empirical_tests = empirical_tests,
                                        empirical_B = empirical_B,
                                        empirical_omnibus_method = empirical_omnibus_method,
                                        empirical_heterogeneity_method = empirical_heterogeneity_method,
                                        empirical_use_as_primary = empirical_use_as_primary,
                                        ci_method = ci_method,
                                        integrated_ci_method = integrated_ci_method,
                                        ci_conf_level = ci_conf_level,
                                        ci_min_boot = ci_min_boot,
                                        orientation = orientation,
                                        pattern_draws = pattern_draws,
                                        analysis_level = analysis_level, seed = seed)),
                   class = "omics_braid_result")
  out
}

#' Run OmicsBraid from externally estimated summary statistics
#'
#' @param effects Data frame with at least `entity`, `omic`, `effect`, and `se`.
#' @param covariance Optional named list of entity-specific covariance matrices.
#' @param omic_order Ordered omics.
#' @param orientation Optional named +1/-1 vector to harmonize omic effect directions.
#' @param equivalence_margin Smallest effect size of interest.
#' @param alpha Significance level.
#' @param state_basis Use local or multiplicity-adjusted inferential states for deterministic braid labels.
#' @param equivalence_p_adjust Multiple-testing method retained for adjusted equivalence/difference evidence.
#' @param trajectory_margin Smallest meaningful standardized-effect change per one-layer transition.
#' @param trend_p_adjust Multiple-testing method retained for adjusted trajectory evidence.
#' @param min_slope Deprecated alias for `trajectory_margin`.
#' @param pattern_draws Monte Carlo draws.
#' @param seed Random seed.
#' @return `omics_braid_result` without sample-level data.
#' @export
run_omics_braid_summary <- function(effects, covariance = NULL, omic_order, orientation = NULL,
                                    equivalence_margin = 0.3, alpha = 0.05,
                                    state_basis = c("local", "adjusted"), equivalence_p_adjust = "BH",
                                    trajectory_margin = 0.15, trend_p_adjust = "BH",
                                    min_slope = NULL, pattern_draws = 2000L, seed = 1L) {
  state_basis <- match.arg(state_basis)
  if (!is.null(min_slope)) {
    trajectory_margin <- min_slope
    .warnf("'min_slope' is deprecated; use 'trajectory_margin'.")
  }
  .validate_effect_table(effects)
  if (anyDuplicated(omic_order)) .stopf("'omic_order' must not contain duplicates.")
  if (!all(unique(effects$omic) %in% omic_order)) .warnf("Some effect-table omics are absent from 'omic_order' and will not contribute to braid classification/plots: %s", paste(setdiff(unique(effects$omic), omic_order), collapse = ", "))
  if (!"p_value" %in% names(effects)) effects$p_value <- 2 * stats::pnorm(abs(effects$effect / effects$se), lower.tail = FALSE)
  if (!"conf_low" %in% names(effects)) effects$conf_low <- effects$effect - stats::qnorm(.975) * effects$se
  if (!"conf_high" %in% names(effects)) effects$conf_high <- effects$effect + stats::qnorm(.975) * effects$se
  if (!"p_adj" %in% names(effects)) effects$p_adj <- .adjust_within_omic(effects$p_value, effects$omic, "BH")
  if (!is.null(orientation)) {
    oriented <- orient_omics(effects, covariance, orientation)
    effects <- oriented$effects; covariance <- oriented$covariance
  }
  integ <- integrate_effects(effects, covariance = covariance)
  eq <- test_equivalence(effects, margin = equivalence_margin, alpha = alpha,
                         p_adjust = equivalence_p_adjust, state_basis = state_basis)
  tr <- test_braid_trend(
    effects, covariance = covariance, omic_order = omic_order,
    trajectory_margin = trajectory_margin, alpha = alpha,
    p_adjust = trend_p_adjust
  )
  cls <- classify_braids(
    eq, omic_order = omic_order, covariance = covariance,
    integrated = integ, trend = tr,
    trajectory_margin = trajectory_margin, alpha = alpha
  )
  probs <- braid_pattern_probabilities(
    effects, covariance = covariance, omic_order = omic_order,
    margin = equivalence_margin, n_draws = pattern_draws,
    trajectory_margin = trajectory_margin, seed = seed
  )
  cls2 <- .summarize_pattern_uncertainty(cls, probs)
  structure(list(data = NULL, effects = effects, covariance = covariance,
                 integrated = integ, equivalence = eq, trend = tr,
                 classification = cls2, pattern_probabilities = probs,
                 settings = list(omic_order = omic_order, orientation = orientation,
                                 equivalence_margin = equivalence_margin,
                                 alpha = alpha, state_basis = state_basis,
                                 equivalence_p_adjust = equivalence_p_adjust,
                                 trajectory_margin = trajectory_margin,
                                 trend_p_adjust = trend_p_adjust,
                                 pattern_draws = pattern_draws,
                                 analysis_level = "summary", seed = seed)),
            class = "omics_braid_result")
}

#' @export
print.omics_braid_result <- function(x, ...) {
  cat("<omics_braid_result>\n")
  cat(" Analysis level:", x$settings$analysis_level, "\n")
  cat(" Entities with integrated estimates:", nrow(x$integrated), "\n")
  cat(" Omic order:", paste(x$settings$omic_order, collapse = " -> "), "\n")
  if (!is.null(x$classification) && nrow(x$classification)) {
    cat(" Pattern counts:\n")
    print(sort(table(x$classification$pattern), decreasing = TRUE))
  }
  invisible(x)
}

#' @export
summary.omics_braid_result <- function(object, ...) {
  list(settings = object$settings,
       n_effects = nrow(object$effects),
       n_integrated = nrow(object$integrated),
       patterns = if (!is.null(object$classification)) sort(table(object$classification$pattern), decreasing = TRUE) else NULL,
       top_evidence = object$integrated[order(object$integrated$p_omnibus_adj), , drop = FALSE][seq_len(min(10L, nrow(object$integrated))), , drop = FALSE],
       top_consensus = object$integrated[order(object$integrated$p_adj), , drop = FALSE][seq_len(min(10L, nrow(object$integrated))), , drop = FALSE])
}
