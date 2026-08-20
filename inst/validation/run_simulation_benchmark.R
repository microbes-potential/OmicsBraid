# OmicsBraid reproducible simulation benchmark
# v0.1.9: inferential calibration + hierarchical braid classification +
# stratified operating characteristics.
#
# Example:
# Rscript run_simulation_benchmark.R --replicates=200 --B=500 --draws=2000 --out=validation_results

validation_lib <- Sys.getenv("OMICSBRAID_VALIDATION_LIB", unset = "")
if (nzchar(validation_lib)) {
  if (!dir.exists(validation_lib)) stop("OMICSBRAID_VALIDATION_LIB does not exist: ", validation_lib, call. = FALSE)
  .libPaths(unique(c(normalizePath(validation_lib, mustWork = TRUE), .libPaths())))
}

suppressPackageStartupMessages(library(OmicsBraid))

expected_version <- Sys.getenv("OMICSBRAID_EXPECTED_VERSION", unset = "")
if (nzchar(expected_version)) {
  installed_version <- as.character(utils::packageVersion("OmicsBraid"))
  if (!identical(installed_version, expected_version)) {
    stop("Validation runner loaded OmicsBraid ", installed_version,
         " but expected ", expected_version, ". Library paths: ",
         paste(.libPaths(), collapse = " | "), call. = FALSE)
  }
}

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", key, "="), "", hit[1])
}
RREP <- as.integer(get_arg("replicates", "100"))
BOOT <- as.integer(get_arg("B", "300"))
DRAWS <- as.integer(get_arg("draws", "1000"))
OUT <- get_arg("out", "OmicsBraid_validation")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

omics <- c("RNA", "Protein", "Metabolite")
truth_patterns <- list(
  joint_null = c(0, 0, 0),
  equal_effect = c(.8, .8, .8),
  concordant_up = c(1.0, .9, .8),
  attenuation = c(1.4, .9, .5),
  amplification = c(.5, .9, 1.4),
  inversion = c(1.0, .8, -1.0),
  buffering = c(1.1, 0, 0),
  emergence = c(0, 0, 1.1)
)
expected_class <- c(
  joint_null = "null_compatible",
  equal_effect = "concordant_increase",
  concordant_up = "concordant_increase",
  attenuation = "attenuation",
  amplification = "amplification",
  inversion = "inversion",
  buffering = "buffering",
  emergence = "emergence"
)

scenarios <- expand.grid(
  n_per_group = c(20L, 40L, 80L),
  rho = c(0, .3, .6),
  missing_rate = c(0, .10),
  residual_distribution = c("normal", "t"),
  replicate = seq_len(RREP),
  stringsAsFactors = FALSE,
  KEEP.OUT.ATTRS = FALSE
)

one_run <- function(z) {
  seed <- 1000000L + as.integer(z$replicate) + 1000L * as.integer(z$n_per_group) +
    10000L * match(z$residual_distribution, c("normal", "t")) +
    as.integer(100 * z$rho) + as.integer(10 * z$missing_rate)
  sim <- simulate_braid_data(
    n_per_group = z$n_per_group,
    omics = omics,
    patterns = truth_patterns,
    rho = z$rho,
    missing_rate = z$missing_rate,
    residual_distribution = z$residual_distribution,
    seed = seed
  )
  fit <- run_omics_braid(
    sim$data, group = "group", reference = "Control", comparison = "Disease",
    omic_order = omics, bootstrap_B = BOOT, bootstrap_shrinkage = .05,
    equivalence_margin = .30, state_basis = "local",
    trajectory_margin = .15, pattern_draws = DRAWS, seed = seed
  )

  integ <- fit$integrated
  cl <- fit$classification
  master <- merge(integ, cl, by = "entity", all.x = TRUE, sort = FALSE)
  master$expected_class <- unname(expected_class[master$entity])

  is_null <- master$entity == "joint_null"
  master$class_correct <- ifelse(
    is_null,
    master$pattern %in% c("null_equivalent", "no_detectable_effect"),
    master$pattern == master$expected_class
  )
  # Hierarchical family correctness distinguishes a safe broader concordant
  # call from a confidently wrong trajectory subtype. For true attenuation or
  # amplification, concordant_increase is a directionally correct but
  # trajectory-unresolved parent class.
  master$family_correct <- master$class_correct
  master$family_correct[master$entity %in% c("attenuation", "amplification")] <-
    master$pattern[master$entity %in% c("attenuation", "amplification")] %in%
    c("attenuation", "amplification", "concordant_increase")
  # Do not treat the opposite trajectory subtype as family-correct.
  master$family_correct[master$entity == "attenuation" & master$pattern == "amplification"] <- FALSE
  master$family_correct[master$entity == "amplification" & master$pattern == "attenuation"] <- FALSE
  master$suggestive_correct <- ifelse(
    is_null,
    master$suggestive_pattern == "null_equivalent",
    master$suggestive_pattern == master$expected_class
  )
  master$recovered_confirmed_or_suggestive <- master$class_correct | master$suggestive_correct
  master$exact_decisive <- master$pattern_status %in% c("confirmed", "no_evidence")
  master$family_decisive <- master$pattern_status %in% c("confirmed", "direction_confirmed", "no_evidence")
  master$wrong_exact_decisive <- master$exact_decisive & !master$class_correct
  master$wrong_family_decisive <- master$family_decisive & !master$family_correct

  # Comparator: what the omnibus/common-effect/Q inference would report if the
  # matched cross-omic dependence were incorrectly ignored.
  indep <- integrate_effects(fit$effects, covariance = NULL)
  ii <- match(master$entity, indep$entity)
  master$p_omnibus_independence <- indep$p_omnibus[ii]
  master$p_consensus_independence <- indep$p_value[ii]
  master$p_heterogeneity_independence <- indep$p_heterogeneity[ii]

  master$n_per_group <- z$n_per_group
  master$rho <- z$rho
  master$missing_rate <- z$missing_rate
  master$residual_distribution <- z$residual_distribution
  master$replicate <- z$replicate

  # Layer-level effect bias and CI coverage against generating shifts.
  ef <- fit$effects
  truth <- sim$truth
  ef <- merge(ef, truth, by = c("entity", "omic"), all.x = TRUE, sort = FALSE)
  eq_keep <- fit$equivalence[, c("entity", "omic", "p_tost", "p_tost_adj",
                                 "equivalent_local", "equivalent_adjusted",
                                 "state_local", "state_adjusted"), drop = FALSE]
  ef <- merge(ef, eq_keep, by = c("entity", "omic"), all.x = TRUE, sort = FALSE)
  ef$bias <- ef$effect - ef$true_effect
  ef$covered95 <- ef$conf_low <= ef$true_effect & ef$conf_high >= ef$true_effect
  ef$n_per_group <- z$n_per_group
  ef$rho <- z$rho
  ef$missing_rate <- z$missing_rate
  ef$residual_distribution <- z$residual_distribution
  ef$replicate <- z$replicate
  list(entity = master, layer = ef)
}

entity_out <- vector("list", nrow(scenarios))
layer_out <- vector("list", nrow(scenarios))
failure_out <- list(); nf <- 0L
for (i in seq_len(nrow(scenarios))) {
  if (i %% 25L == 0L) message("Scenario run ", i, " / ", nrow(scenarios))
  ans <- try(one_run(scenarios[i, , drop = FALSE]), silent = TRUE)
  if (inherits(ans, "try-error")) {
    nf <- nf + 1L
    failure_out[[nf]] <- data.frame(
      scenario_row = i,
      scenarios[i, , drop = FALSE],
      error = as.character(ans),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    next
  }
  entity_out[[i]] <- ans$entity
  layer_out[[i]] <- ans$layer
}

entity_ok <- entity_out[vapply(entity_out, function(x) is.data.frame(x) && "entity" %in% names(x), logical(1))]
layer_ok <- layer_out[vapply(layer_out, function(x) is.data.frame(x) && "entity" %in% names(x), logical(1))]
if (!length(entity_ok) || !length(layer_ok)) stop("No simulation scenarios completed successfully.", call. = FALSE)
entity_results <- do.call(rbind, entity_ok)
layer_results <- do.call(rbind, layer_ok)
failures <- if (length(failure_out)) do.call(rbind, failure_out) else data.frame(
  scenario_row = integer(), n_per_group = integer(), rho = numeric(),
  missing_rate = numeric(), residual_distribution = character(), replicate = integer(),
  error = character(), stringsAsFactors = FALSE
)
utils::write.csv(entity_results, file.path(OUT, "entity_level_results.csv"), row.names = FALSE)
utils::write.csv(layer_results, file.path(OUT, "layer_level_results.csv"), row.names = FALSE)
utils::write.csv(failures, file.path(OUT, "scenario_failures.csv"), row.names = FALSE)

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

wilson <- function(x, conf = .95) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (!n) return(c(value = NA_real_, n = 0, mc_se = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  k <- sum(x)
  p <- k / n
  z <- stats::qnorm(1 - (1 - conf) / 2)
  z2 <- z^2
  den <- 1 + z2 / n
  center <- (p + z2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z2 / (4 * n^2)) / den
  c(value = p, n = n, mc_se = sqrt(p * (1 - p) / n),
    ci_low = max(0, center - half), ci_high = min(1, center + half))
}

metric_row <- function(name, x) {
  z <- wilson(x)
  data.frame(metric = name, value = unname(z["value"]), n_evaluable = as.integer(z["n"]),
             mc_se = unname(z["mc_se"]), mc_ci_low = unname(z["ci_low"]),
             mc_ci_high = unname(z["ci_high"]), stringsAsFactors = FALSE)
}

make_core_metrics <- function(er, lr) {
  null <- er$entity == "joint_null"
  equal <- er$entity == "equal_effect"
  inv <- er$entity == "inversion"
  att <- er$entity == "attenuation"
  amp <- er$entity == "amplification"
  rows <- list(
    metric_row("omnibus_type_I_at_0.05_joint_null", er$p_omnibus[null] < .05),
    metric_row("consensus_type_I_at_0.05_joint_null", er$p_value[null] < .05),
    metric_row("Q_type_I_at_0.05_equal_effect", er$p_heterogeneity[equal] < .05),
    metric_row("omnibus_power_inversion", er$p_omnibus[inv] < .05),
    metric_row("Q_power_inversion", er$p_heterogeneity[inv] < .05),
    metric_row("trend_false_meaningful_call_equal_effect", er$trend_state_local[equal] %in% c("attenuation", "amplification")),
    metric_row("trend_flat_confirmation_equal_effect", er$trend_state_local[equal] == "flat"),
    metric_row("trend_power_attenuation", er$trend_state_local[att] == "attenuation"),
    metric_row("trend_power_amplification", er$trend_state_local[amp] == "amplification"),
    metric_row("strict_overall_pattern_accuracy", er$class_correct),
    metric_row("hierarchical_family_accuracy", er$family_correct),
    metric_row("exact_decisive_call_rate", er$exact_decisive),
    metric_row("family_decisive_call_rate", er$family_decisive),
    metric_row("wrong_exact_decisive_call_rate", er$wrong_exact_decisive),
    metric_row("wrong_family_decisive_call_rate", er$wrong_family_decisive),
    metric_row("confirmed_or_suggestive_recovery_rate", er$recovered_confirmed_or_suggestive),
    metric_row("joint_null_compatible_primary_rate", er$class_correct[null]),
    metric_row("joint_null_equivalence_confirmation_rate", er$pattern[null] == "null_equivalent"),
    metric_row("joint_null_no_detectable_effect_rate", er$pattern[null] == "no_detectable_effect"),
    metric_row("layer_95CI_coverage", lr$covered95),
    metric_row("true_zero_layer_equivalence_confirmation_rate", lr$equivalent_local[abs(lr$true_effect) < 1e-12]),
    metric_row("false_equivalence_rate_for_true_non_negligible_layers", lr$equivalent_local[abs(lr$true_effect) > .30 + 1e-12]),
    metric_row("equal_effect_consensus_95CI_coverage", er$conf_low[equal] <= .8 & er$conf_high[equal] >= .8),
    metric_row("omnibus_type_I_joint_null_if_independence_assumed", er$p_omnibus_independence[null] < .05),
    metric_row("Q_type_I_equal_effect_if_independence_assumed", er$p_heterogeneity_independence[equal] < .05)
  )
  do.call(rbind, rows)
}

summary_metrics <- make_core_metrics(entity_results, layer_results)
# Scenario execution itself is also part of operating reliability.
summary_metrics <- rbind(
  summary_metrics,
  data.frame(metric = "scenario_failure_rate", value = nrow(failures) / nrow(scenarios),
             n_evaluable = nrow(scenarios), mc_se = NA_real_, mc_ci_low = NA_real_, mc_ci_high = NA_real_,
             stringsAsFactors = FALSE)
)
utils::write.csv(summary_metrics, file.path(OUT, "core_validation_metrics.csv"), row.names = FALSE)

make_pattern_metrics <- function(d) {
  exact <- d$exact_decisive %in% TRUE
  fam <- d$family_decisive %in% TRUE
  data.frame(
    expected_class = d$expected_class[1],
    n = nrow(d),
    strict_accuracy = mean(d$class_correct, na.rm = TRUE),
    hierarchical_family_accuracy = mean(d$family_correct, na.rm = TRUE),
    exact_decisive_rate = mean(d$exact_decisive, na.rm = TRUE),
    family_decisive_rate = mean(d$family_decisive, na.rm = TRUE),
    accuracy_if_exact_decisive = if (any(exact, na.rm = TRUE)) mean(d$class_correct[exact], na.rm = TRUE) else NA_real_,
    family_accuracy_if_family_decisive = if (any(fam, na.rm = TRUE)) mean(d$family_correct[fam], na.rm = TRUE) else NA_real_,
    wrong_exact_decisive_rate = mean(d$wrong_exact_decisive, na.rm = TRUE),
    wrong_family_decisive_rate = mean(d$wrong_family_decisive, na.rm = TRUE),
    suggestive_accuracy = mean(d$suggestive_correct, na.rm = TRUE),
    confirmed_or_suggestive_recovery = mean(d$recovered_confirmed_or_suggestive, na.rm = TRUE),
    confirmed_status_rate = mean(d$pattern_status == "confirmed", na.rm = TRUE),
    direction_confirmed_status_rate = mean(d$pattern_status == "direction_confirmed", na.rm = TRUE),
    no_evidence_status_rate = mean(d$pattern_status == "no_evidence", na.rm = TRUE),
    unresolved_status_rate = mean(d$pattern_status %in% c("unresolved", "insufficient"), na.rm = TRUE),
    mean_deterministic_support = safe_mean(d$deterministic_pattern_support),
    mean_suggestive_support = safe_mean(d$suggestive_pattern_support),
    mean_pattern_entropy = safe_mean(d$pattern_entropy),
    stringsAsFactors = FALSE
  )
}
pattern_metrics <- do.call(rbind, lapply(split(entity_results, entity_results$expected_class), make_pattern_metrics))
utils::write.csv(pattern_metrics, file.path(OUT, "pattern_operating_characteristics.csv"), row.names = FALSE)

# Trend-specific operating characteristics.
trend_metrics <- data.frame(
  truth = c("equal_effect", "concordant_up", "attenuation", "amplification"),
  expected_trend = c("flat", "practically_flat_or_unresolved", "attenuation", "amplification"),
  n = c(sum(entity_results$entity == "equal_effect"), sum(entity_results$entity == "concordant_up"),
        sum(entity_results$entity == "attenuation"), sum(entity_results$entity == "amplification")),
  correct_or_acceptable = c(
    mean(entity_results$trend_state_local[entity_results$entity == "equal_effect"] == "flat", na.rm = TRUE),
    mean(entity_results$trend_state_local[entity_results$entity == "concordant_up"] %in% c("flat", "uncertain"), na.rm = TRUE),
    mean(entity_results$trend_state_local[entity_results$entity == "attenuation"] == "attenuation", na.rm = TRUE),
    mean(entity_results$trend_state_local[entity_results$entity == "amplification"] == "amplification", na.rm = TRUE)
  ),
  wrong_opposite_trend = c(
    mean(entity_results$trend_state_local[entity_results$entity == "equal_effect"] %in% c("attenuation", "amplification"), na.rm = TRUE),
    mean(entity_results$trend_state_local[entity_results$entity == "concordant_up"] == "amplification", na.rm = TRUE),
    mean(entity_results$trend_state_local[entity_results$entity == "attenuation"] == "amplification", na.rm = TRUE),
    mean(entity_results$trend_state_local[entity_results$entity == "amplification"] == "attenuation", na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(trend_metrics, file.path(OUT, "trend_operating_characteristics.csv"), row.names = FALSE)

# Confusion matrix for primary inferential labels. The true zero-generating
# scenario is deliberately called null_compatible because no-detectable-effect
# and demonstrated-equivalence are distinct but both scientifically valid.
conf <- with(entity_results, table(expected = expected_class, observed = pattern, useNA = "ifany"))
utils::write.csv(as.data.frame.matrix(conf), file.path(OUT, "braid_confusion_matrix.csv"))

# Generic grouped summaries for the factors that define the benchmark.
group_pattern_summary <- function(er, vars) {
  key <- interaction(er[, vars, drop = FALSE], drop = TRUE, lex.order = TRUE)
  pieces <- lapply(split(seq_len(nrow(er)), key), function(ii) {
    d <- er[ii, , drop = FALSE]
    # make_pattern_metrics() already contains expected_class; avoid duplicate
    # column names in grouped outputs.
    extra_vars <- setdiff(vars, "expected_class")
    base <- if (length(extra_vars)) d[1, extra_vars, drop = FALSE] else data.frame()
    cbind(base, make_pattern_metrics(d))
  })
  do.call(rbind, pieces)
}
utils::write.csv(group_pattern_summary(entity_results, c("expected_class", "n_per_group")),
                 file.path(OUT, "pattern_metrics_by_sample_size.csv"), row.names = FALSE)
utils::write.csv(group_pattern_summary(entity_results, c("expected_class", "rho")),
                 file.path(OUT, "pattern_metrics_by_correlation.csv"), row.names = FALSE)
utils::write.csv(group_pattern_summary(entity_results, c("expected_class", "missing_rate")),
                 file.path(OUT, "pattern_metrics_by_missingness.csv"), row.names = FALSE)
utils::write.csv(group_pattern_summary(entity_results, c("expected_class", "residual_distribution")),
                 file.path(OUT, "pattern_metrics_by_distribution.csv"), row.names = FALSE)
utils::write.csv(group_pattern_summary(entity_results,
                                       c("expected_class", "n_per_group", "rho", "missing_rate", "residual_distribution")),
                 file.path(OUT, "scenario_stratified_pattern_metrics.csv"), row.names = FALSE)

core_by_factor <- function(er, lr, var) {
  vals <- unique(er[[var]])
  pieces <- lapply(vals, function(v) {
    er0 <- er[er[[var]] == v, , drop = FALSE]
    lr0 <- lr[lr[[var]] == v, , drop = FALSE]
    z <- make_core_metrics(er0, lr0)
    z$factor <- var
    z$level <- as.character(v)
    z
  })
  z <- do.call(rbind, pieces)
  z[, c("factor", "level", "metric", "value", "n_evaluable", "mc_se", "mc_ci_low", "mc_ci_high")]
}
utils::write.csv(core_by_factor(entity_results, layer_results, "n_per_group"),
                 file.path(OUT, "core_metrics_by_sample_size.csv"), row.names = FALSE)
utils::write.csv(core_by_factor(entity_results, layer_results, "rho"),
                 file.path(OUT, "core_metrics_by_correlation.csv"), row.names = FALSE)
utils::write.csv(core_by_factor(entity_results, layer_results, "missing_rate"),
                 file.path(OUT, "core_metrics_by_missingness.csv"), row.names = FALSE)
utils::write.csv(core_by_factor(entity_results, layer_results, "residual_distribution"),
                 file.path(OUT, "core_metrics_by_distribution.csv"), row.names = FALSE)

# Preserve everything needed for audit/re-analysis.
saveRDS(list(
  settings = list(replicates = RREP, bootstrap_B = BOOT, pattern_draws = DRAWS,
                  equivalence_margin = .30, trajectory_margin = .15),
  scenarios = scenarios,
  failures = failures,
  entity_results = entity_results,
  layer_results = layer_results,
  core_metrics = summary_metrics,
  pattern_metrics = pattern_metrics,
  trend_metrics = trend_metrics
), file.path(OUT, "validation_bundle.rds"))

capture.output(sessionInfo(), file = file.path(OUT, "sessionInfo.txt"))
message("Validation benchmark written to: ", normalizePath(OUT, mustWork = FALSE))
