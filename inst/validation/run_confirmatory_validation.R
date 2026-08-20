# OmicsBraid v0.2.1 targeted confirmatory validation (I/O-resilience patch; statistical design unchanged)
#
# Purpose:
#   1) precisely calibrate the omnibus and Q_omics Type-I error by sample size;
#   2) compare analytic and subject-bootstrap confidence-interval coverage;
#   3) quantify heavy-tail robustness;
#   4) estimate attenuation/amplification trend power curves;
#   5) estimate equivalence power for true-zero buffering/emergence layers;
#   6) compare matched covariance against an incorrect independence assumption;
#   7) audit decisive vs unresolved braid classification.
#
# Example:
# Rscript run_confirmatory_validation.R --replicates=500 --B=500 --bca-replicates=50 --cores=2 --out=confirmatory_results

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
    stop("Loaded OmicsBraid ", installed_version, " but expected ", expected_version, call. = FALSE)
  }
}

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", key, "="), "", hit[1])
}
RREP <- as.integer(get_arg("replicates", "500"))
BOOT <- as.integer(get_arg("B", "500"))
BCA_REP <- as.integer(get_arg("bca-replicates", "50"))
CORES <- as.integer(get_arg("cores", "2"))
CHUNK <- as.integer(get_arg("chunk", "20"))
OUT <- get_arg("out", "OmicsBraid_confirmatory_validation")
if (!is.finite(RREP) || RREP < 50) stop("--replicates must be >= 50", call. = FALSE)
if (!is.finite(BOOT) || BOOT < 100) stop("--B must be >= 100", call. = FALSE)
if (!is.finite(CORES) || CORES < 1) CORES <- 1L
if (!is.finite(CHUNK) || CHUNK < 1) CHUNK <- 20L
if (!is.finite(BCA_REP) || BCA_REP < 0) BCA_REP <- 0L
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
checkpoint_dir <- file.path(OUT, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

checkpoint_is_valid <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0) return(FALSE)
  z <- try(readRDS(path), silent = TRUE)
  if (inherits(z, "try-error")) return(FALSE)
  is.list(z) && length(z) > 0L
}

# Write a checkpoint to a temporary file on the SAME filesystem and rename it
# only after serialization succeeds. This prevents partial files from being
# mistaken for completed chunks after an interrupted write.
atomic_save_rds <- function(object, path, attempts = 5L) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  for (attempt in seq_len(attempts)) {
    tmp <- tempfile(pattern = paste0(".", basename(path), ".tmp_"),
                    tmpdir = dirname(path), fileext = ".rds")
    ok <- tryCatch({
      saveRDS(object, tmp, compress = FALSE)
      if (!file.exists(tmp) || file.info(tmp)$size <= 0) stop("temporary checkpoint was not written")
      if (file.exists(path)) unlink(path, force = TRUE)
      if (!file.rename(tmp, path)) stop("atomic checkpoint rename failed")
      checkpoint_is_valid(path)
    }, error = function(e) FALSE)
    unlink(tmp, force = TRUE)
    if (isTRUE(ok)) return(invisible(TRUE))
    Sys.sleep(0.25 * attempt)
  }
  stop("Unable to write a valid checkpoint after ", attempts, " attempts: ", path, call. = FALSE)
}

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

# Main confirmatory grid deliberately focuses on sample size and tail behavior.
# Correlation and modest missingness were already explored in v0.1.9; here rho
# is fixed at a realistic moderate dependence so each n x distribution cell has
# enough Monte Carlo precision for calibration.
scenarios <- expand.grid(
  n_per_group = c(20L, 40L, 80L, 160L, 320L),
  residual_distribution = c("normal", "t"),
  replicate = seq_len(RREP),
  stringsAsFactors = FALSE,
  KEEP.OUT.ATTRS = FALSE
)
scenarios$rho <- 0.30
scenarios$missing_rate <- 0
scenarios$scenario_id <- seq_len(nrow(scenarios))

seed_for <- function(z) {
  20000000L + as.integer(z$replicate) +
    10000L * as.integer(z$n_per_group) +
    1000000L * match(z$residual_distribution, c("normal", "t"))
}

one_run <- function(z) {
  seed <- seed_for(z)
  sim <- simulate_braid_data(
    n_per_group = z$n_per_group,
    omics = omics,
    patterns = truth_patterns,
    rho = z$rho,
    missing_rate = z$missing_rate,
    residual_distribution = z$residual_distribution,
    seed = seed
  )
  eff <- estimate_effects(sim$data, "group", "Control", "Disease", conf_level = .95)
  cv <- bootstrap_effect_covariance(
    sim$data, "group", "Control", "Disease",
    effects = eff, B = BOOT, seed = seed,
    shrinkage = .05, min_complete = max(50L, floor(BOOT * .50)),
    stratified = TRUE
  )
  integ <- integrate_effects(eff, covariance = cv, conf_level = .95)
  eq <- test_equivalence(eff, margin = .30, alpha = .05, state_basis = "local")
  tr <- test_braid_trend(
    eff, covariance = cv, omic_order = omics,
    trajectory_margin = .15, alpha = .05
  )
  cl <- classify_braids(
    eq, omic_order = omics, covariance = cv,
    integrated = integ, trend = tr,
    trajectory_margin = .15, alpha = .05
  )

  # Robust interval comparisons reuse exactly the same matched-subject bootstrap
  # that estimated cross-omic dependence; inference p-values remain analytic.
  ci_pct <- bootstrap_effect_intervals(
    eff, cv, method = "percentile", conf_level = .95,
    min_boot = max(50L, floor(BOOT * .50))
  )
  ci_basic <- bootstrap_effect_intervals(
    eff, cv, method = "basic", conf_level = .95,
    min_boot = max(50L, floor(BOOT * .50))
  )
  cons_pct <- bootstrap_consensus_intervals(
    eff, cv, integrated = integ, method = "percentile",
    conf_level = .95, min_boot = max(50L, floor(BOOT * .50))
  )

  # BCa is substantially more expensive because it needs leave-one-subject-out
  # acceleration. Validate it on a prespecified subset of heavy-tail scenarios.
  ci_bca <- NULL
  if (BCA_REP > 0L && z$residual_distribution == "t" &&
      z$n_per_group %in% c(40L, 80L, 160L) && z$replicate <= BCA_REP) {
    ci_bca <- bootstrap_effect_intervals(
      eff, cv, method = "bca", conf_level = .95,
      min_boot = max(50L, floor(BOOT * .50)),
      data = sim$data, group = "group",
      reference = "Control", comparison = "Disease"
    )
  }

  indep <- integrate_effects(eff, covariance = NULL, conf_level = .95)

  # classify_braids() already carries the trend diagnostics needed for
  # operating-characteristic summaries; avoid duplicate .x/.y columns.
  master <- merge(integ, cl, by = "entity", all.x = TRUE, sort = FALSE)
  master$expected_class <- unname(expected_class[master$entity])
  is_null <- master$entity == "joint_null"
  master$class_correct <- ifelse(
    is_null,
    master$pattern %in% c("null_equivalent", "no_detectable_effect"),
    master$pattern == master$expected_class
  )
  master$family_correct <- master$class_correct
  ktraj <- master$entity %in% c("attenuation", "amplification")
  master$family_correct[ktraj] <- master$pattern[ktraj] %in%
    c("attenuation", "amplification", "concordant_increase")
  master$family_correct[master$entity == "attenuation" & master$pattern == "amplification"] <- FALSE
  master$family_correct[master$entity == "amplification" & master$pattern == "attenuation"] <- FALSE
  master$exact_decisive <- master$pattern_status %in% c("confirmed", "no_evidence")
  master$family_decisive <- master$pattern_status %in% c("confirmed", "direction_confirmed", "no_evidence")
  master$wrong_exact_decisive <- master$exact_decisive & !master$class_correct
  master$wrong_family_decisive <- master$family_decisive & !master$family_correct

  ii <- match(master$entity, indep$entity)
  master$p_omnibus_independence <- indep$p_omnibus[ii]
  master$p_consensus_independence <- indep$p_value[ii]
  master$p_heterogeneity_independence <- indep$p_heterogeneity[ii]
  ip <- match(master$entity, cons_pct$entity)
  master$consensus_pct_low <- cons_pct$conf_low_boot[ip]
  master$consensus_pct_high <- cons_pct$conf_high_boot[ip]
  master$n_per_group <- z$n_per_group
  master$residual_distribution <- z$residual_distribution
  master$rho <- z$rho
  master$replicate <- z$replicate
  master$scenario_id <- z$scenario_id

  ef <- merge(eff, sim$truth, by = c("entity", "omic"), all.x = TRUE, sort = FALSE)
  eqk <- eq[, c("entity", "omic", "equivalent_local", "state_local"), drop = FALSE]
  ef <- merge(ef, eqk, by = c("entity", "omic"), all.x = TRUE, sort = FALSE)
  pkey <- paste(ef$entity, ef$omic, sep = "\r")
  ikey <- paste(ci_pct$entity, ci_pct$omic, sep = "\r")
  ibasic <- paste(ci_basic$entity, ci_basic$omic, sep = "\r")
  jj <- match(pkey, ikey)
  jb <- match(pkey, ibasic)
  ef$analytic_covered95 <- ef$conf_low <= ef$true_effect & ef$conf_high >= ef$true_effect
  ef$percentile_low <- ci_pct$conf_low_boot[jj]
  ef$percentile_high <- ci_pct$conf_high_boot[jj]
  ef$percentile_covered95 <- ef$percentile_low <= ef$true_effect & ef$percentile_high >= ef$true_effect
  ef$basic_low <- ci_basic$conf_low_boot[jb]
  ef$basic_high <- ci_basic$conf_high_boot[jb]
  ef$basic_covered95 <- ef$basic_low <= ef$true_effect & ef$basic_high >= ef$true_effect
  ef$bca_low <- NA_real_; ef$bca_high <- NA_real_; ef$bca_covered95 <- NA
  if (!is.null(ci_bca) && nrow(ci_bca)) {
    kb <- paste(ci_bca$entity, ci_bca$omic, sep = "\r")
    jbc <- match(pkey, kb)
    ef$bca_low <- ci_bca$conf_low_boot[jbc]
    ef$bca_high <- ci_bca$conf_high_boot[jbc]
    ef$bca_covered95 <- ef$bca_low <= ef$true_effect & ef$bca_high >= ef$true_effect
  }
  ef$n_per_group <- z$n_per_group
  ef$residual_distribution <- z$residual_distribution
  ef$rho <- z$rho
  ef$replicate <- z$replicate
  ef$scenario_id <- z$scenario_id

  list(entity = master, layer = ef)
}

# Chunked checkpointing makes the long confirmatory run safely resumable.
ids <- scenarios$scenario_id
chunks <- split(ids, ceiling(seq_along(ids) / CHUNK))
for (cc in seq_along(chunks)) {
  f <- file.path(checkpoint_dir, sprintf("chunk_%04d.rds", cc))
  if (checkpoint_is_valid(f)) next
  if (file.exists(f)) {
    warning("Removing unreadable/incomplete checkpoint before recomputation: ", f, call. = FALSE)
    unlink(f, force = TRUE)
  }
  idx <- chunks[[cc]]
  run_one_id <- function(i) {
    z <- scenarios[scenarios$scenario_id == i, , drop = FALSE]
    ans <- try(one_run(z), silent = TRUE)
    if (inherits(ans, "try-error")) {
      return(list(ok = FALSE, scenario = z, error = as.character(ans)))
    }
    list(ok = TRUE, scenario = z, result = ans)
  }
  if (.Platform$OS.type != "windows" && CORES > 1L) {
    ans <- parallel::mclapply(idx, run_one_id, mc.cores = CORES, mc.preschedule = FALSE)
  } else {
    ans <- lapply(idx, run_one_id)
  }
  atomic_save_rds(ans, f)
  message("Completed confirmatory chunk ", cc, " / ", length(chunks),
          " (scenario rows through ", max(idx), " / ", nrow(scenarios), ")")
}

chunk_files <- sort(list.files(checkpoint_dir, pattern = "^chunk_[0-9]+\\.rds$", full.names = TRUE))
all_runs <- unlist(lapply(chunk_files, readRDS), recursive = FALSE)
oks <- vapply(all_runs, function(x) isTRUE(x$ok), logical(1))
fail_runs <- all_runs[!oks]
good_runs <- all_runs[oks]
if (!length(good_runs)) stop("No confirmatory scenarios completed successfully.", call. = FALSE)
entity_results <- do.call(rbind, lapply(good_runs, function(x) x$result$entity))
layer_results <- do.call(rbind, lapply(good_runs, function(x) x$result$layer))
failures <- if (length(fail_runs)) do.call(rbind, lapply(fail_runs, function(x) {
  cbind(x$scenario, error = x$error, stringsAsFactors = FALSE)
})) else data.frame(
  n_per_group = integer(), residual_distribution = character(), replicate = integer(),
  rho = numeric(), missing_rate = numeric(), scenario_id = integer(), error = character()
)
utils::write.csv(entity_results, file.path(OUT, "entity_level_confirmatory.csv"), row.names = FALSE)
utils::write.csv(layer_results, file.path(OUT, "layer_level_confirmatory.csv"), row.names = FALSE)
utils::write.csv(failures, file.path(OUT, "scenario_failures.csv"), row.names = FALSE)

wilson <- function(x, conf = .95) {
  x <- x[!is.na(x)]; n <- length(x)
  if (!n) return(c(value = NA_real_, n = 0, mc_se = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  k <- sum(x); p <- k / n; z <- stats::qnorm(1 - (1 - conf) / 2); z2 <- z^2
  den <- 1 + z2 / n
  center <- (p + z2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z2 / (4 * n^2)) / den
  c(value = p, n = n, mc_se = sqrt(p * (1 - p) / n),
    ci_low = max(0, center - half), ci_high = min(1, center + half))
}
metric <- function(name, x, target = NA_real_) {
  z <- wilson(x)
  data.frame(metric = name, value = unname(z["value"]), n_evaluable = as.integer(z["n"]),
             mc_se = unname(z["mc_se"]), mc_ci_low = unname(z["ci_low"]),
             mc_ci_high = unname(z["ci_high"]), target = target,
             target_inside_mc_ci = if (is.finite(target)) unname(z["ci_low"] <= target && z["ci_high"] >= target) else NA,
             stringsAsFactors = FALSE)
}

summarize_cell <- function(er, lr, n, dist) {
  e <- er[er$n_per_group == n & er$residual_distribution == dist, , drop = FALSE]
  l <- lr[lr$n_per_group == n & lr$residual_distribution == dist, , drop = FALSE]
  null <- e$entity == "joint_null"; equal <- e$entity == "equal_effect"
  att <- e$entity == "attenuation"; amp <- e$entity == "amplification"
  inv <- e$entity == "inversion"
  z0 <- abs(l$true_effect) < 1e-12
  nonneg <- abs(l$true_effect) > .30 + 1e-12
  equal_cov_pct <- e$consensus_pct_low[equal] <= .8 & e$consensus_pct_high[equal] >= .8
  rows <- rbind(
    metric("omnibus_type_I", e$p_omnibus[null] < .05, .05),
    metric("consensus_type_I", e$p_value[null] < .05, .05),
    metric("Q_type_I_equal_effect", e$p_heterogeneity[equal] < .05, .05),
    metric("omnibus_power_inversion", e$p_omnibus[inv] < .05),
    metric("Q_power_inversion", e$p_heterogeneity[inv] < .05),
    metric("trend_false_meaningful_equal", e$trend_state_local[equal] %in% c("attenuation", "amplification"), .05),
    metric("trend_power_attenuation", e$trend_state_local[att] == "attenuation"),
    metric("trend_power_amplification", e$trend_state_local[amp] == "amplification"),
    metric("wrong_exact_decisive", e$wrong_exact_decisive, .05),
    metric("wrong_family_decisive", e$wrong_family_decisive, .05),
    metric("analytic_layer_95CI_coverage", l$analytic_covered95, .95),
    metric("percentile_layer_95CI_coverage", l$percentile_covered95, .95),
    metric("basic_layer_95CI_coverage", l$basic_covered95, .95),
    metric("bca_layer_95CI_coverage", l$bca_covered95, .95),
    metric("true_zero_equivalence_power", l$equivalent_local[z0]),
    metric("false_equivalence_non_negligible", l$equivalent_local[nonneg], 0),
    metric("analytic_equal_consensus_95CI_coverage", e$conf_low[equal] <= .8 & e$conf_high[equal] >= .8, .95),
    metric("percentile_equal_consensus_95CI_coverage", equal_cov_pct, .95),
    metric("omnibus_type_I_if_independence_assumed", e$p_omnibus_independence[null] < .05, .05),
    metric("Q_type_I_if_independence_assumed", e$p_heterogeneity_independence[equal] < .05, .05)
  )
  rows$n_per_group <- n; rows$residual_distribution <- dist
  rows[, c("n_per_group", "residual_distribution", "metric", "value", "n_evaluable", "mc_se", "mc_ci_low", "mc_ci_high", "target", "target_inside_mc_ci")]
}

cells <- do.call(rbind, lapply(c(20L,40L,80L,160L,320L), function(n) {
  do.call(rbind, lapply(c("normal","t"), function(d) summarize_cell(entity_results, layer_results, n, d)))
}))
utils::write.csv(cells, file.path(OUT, "confirmatory_metrics_by_sample_size_distribution.csv"), row.names = FALSE)

# More readable focused tables.
pick_metrics <- function(names) cells[cells$metric %in% names, , drop = FALSE]
utils::write.csv(pick_metrics(c("omnibus_type_I", "Q_type_I_equal_effect", "omnibus_type_I_if_independence_assumed", "Q_type_I_if_independence_assumed")),
                 file.path(OUT, "type1_calibration_by_sample_size.csv"), row.names = FALSE)
utils::write.csv(pick_metrics(c("analytic_layer_95CI_coverage", "percentile_layer_95CI_coverage", "basic_layer_95CI_coverage", "bca_layer_95CI_coverage", "analytic_equal_consensus_95CI_coverage", "percentile_equal_consensus_95CI_coverage")),
                 file.path(OUT, "ci_coverage_by_sample_size.csv"), row.names = FALSE)
utils::write.csv(pick_metrics(c("trend_false_meaningful_equal", "trend_power_attenuation", "trend_power_amplification")),
                 file.path(OUT, "trend_power_by_sample_size.csv"), row.names = FALSE)
utils::write.csv(pick_metrics(c("true_zero_equivalence_power", "false_equivalence_non_negligible")),
                 file.path(OUT, "equivalence_power_by_sample_size.csv"), row.names = FALSE)
utils::write.csv(pick_metrics(c("wrong_exact_decisive", "wrong_family_decisive")),
                 file.path(OUT, "classification_safety_by_sample_size.csv"), row.names = FALSE)

# Aggregate heavy-tail CI comparison across n.
ci_dist <- do.call(rbind, lapply(c("normal", "t"), function(d) {
  l <- layer_results[layer_results$residual_distribution == d, , drop = FALSE]
  rbind(
    cbind(residual_distribution = d, metric("analytic_layer_95CI_coverage", l$analytic_covered95, .95)),
    cbind(residual_distribution = d, metric("percentile_layer_95CI_coverage", l$percentile_covered95, .95)),
    cbind(residual_distribution = d, metric("basic_layer_95CI_coverage", l$basic_covered95, .95)),
    cbind(residual_distribution = d, metric("bca_layer_95CI_coverage", l$bca_covered95, .95))
  )
}))
utils::write.csv(ci_dist, file.path(OUT, "ci_coverage_by_distribution.csv"), row.names = FALSE)

# Acceptance-oriented summary: calibration criteria are evaluated against their
# Monte Carlo intervals, not by requiring point estimates to equal nominal values.
acceptance <- rbind(
  data.frame(
    criterion = "Omnibus Type-I calibration",
    result = if (all(cells$target_inside_mc_ci[cells$metric == "omnibus_type_I"] %in% TRUE)) "PASS" else "REVIEW",
    rule = "Nominal 0.05 lies inside each n x distribution Monte Carlo 95% interval",
    stringsAsFactors = FALSE
  ),
  data.frame(
    criterion = "Q_omics Type-I calibration",
    result = if (all(cells$target_inside_mc_ci[cells$metric == "Q_type_I_equal_effect"] %in% TRUE)) "PASS" else "REVIEW",
    rule = "Nominal 0.05 lies inside each n x distribution Monte Carlo 95% interval",
    stringsAsFactors = FALSE
  ),
  data.frame(
    criterion = "Bootstrap percentile layer-CI coverage",
    result = if (all(cells$target_inside_mc_ci[cells$metric == "percentile_layer_95CI_coverage"] %in% TRUE)) "PASS" else "REVIEW",
    rule = "Nominal 0.95 lies inside each n x distribution Monte Carlo 95% interval",
    stringsAsFactors = FALSE
  ),
  data.frame(
    criterion = "Wrong family-decisive classification",
    result = if (max(cells$value[cells$metric == "wrong_family_decisive"], na.rm = TRUE) <= .05) "PASS" else "REVIEW",
    rule = "Point estimate <= 0.05 in every n x distribution cell",
    stringsAsFactors = FALSE
  ),
  data.frame(
    criterion = "False equivalence of non-negligible effects",
    result = if (max(cells$value[cells$metric == "false_equivalence_non_negligible"], na.rm = TRUE) <= .01) "PASS" else "REVIEW",
    rule = "Point estimate <= 0.01 in every n x distribution cell",
    stringsAsFactors = FALSE
  ),
  data.frame(
    criterion = "Scenario execution reliability",
    result = if (nrow(failures) == 0L) "PASS" else "REVIEW",
    rule = "No failed scenario rows",
    stringsAsFactors = FALSE
  )
)
utils::write.csv(acceptance, file.path(OUT, "acceptance_summary.csv"), row.names = FALSE)

# Publication-style validation plots.
if (requireNamespace("ggplot2", quietly = TRUE)) {
  make_plot <- function(dat, metric_name, ylab, target = NULL, file) {
    d <- dat[dat$metric == metric_name & is.finite(dat$value), , drop = FALSE]
    p <- ggplot2::ggplot(d, ggplot2::aes(x = n_per_group, y = value,
                                         linetype = residual_distribution,
                                         group = residual_distribution)) +
      ggplot2::geom_line() + ggplot2::geom_point() +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = mc_ci_low, ymax = mc_ci_high), width = 0) +
      ggplot2::scale_x_continuous(breaks = c(20,40,80,160,320)) +
      ggplot2::labs(x = "Sample size per group", y = ylab,
                    linetype = "Residual distribution") +
      ggplot2::theme_classic()
    if (!is.null(target)) p <- p + ggplot2::geom_hline(yintercept = target, linetype = 3)
    ggplot2::ggsave(file.path(OUT, file), p, width = 7, height = 5)
  }
  make_plot(cells, "omnibus_type_I", "Omnibus Type-I error", .05, "Fig_validation_omnibus_type1.pdf")
  make_plot(cells, "Q_type_I_equal_effect", "Qomics Type-I error", .05, "Fig_validation_Q_type1.pdf")
  make_plot(cells, "trend_power_attenuation", "Attenuation trend power", NULL, "Fig_validation_attenuation_power.pdf")
  make_plot(cells, "trend_power_amplification", "Amplification trend power", NULL, "Fig_validation_amplification_power.pdf")
  make_plot(cells, "true_zero_equivalence_power", "Equivalence power for true-zero layers", NULL, "Fig_validation_equivalence_power.pdf")

  dci <- cells[cells$metric %in% c("analytic_layer_95CI_coverage", "percentile_layer_95CI_coverage") & is.finite(cells$value), , drop = FALSE]
  dci$method <- ifelse(dci$metric == "analytic_layer_95CI_coverage", "Analytic", "Percentile bootstrap")
  p <- ggplot2::ggplot(dci, ggplot2::aes(x = n_per_group, y = value, linetype = method,
                                         shape = residual_distribution,
                                         group = interaction(method, residual_distribution))) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = .95, linetype = 3) +
    ggplot2::scale_x_continuous(breaks = c(20,40,80,160,320)) +
    ggplot2::labs(x = "Sample size per group", y = "95% CI coverage",
                  linetype = "Interval", shape = "Residual distribution") +
    ggplot2::theme_classic()
  ggplot2::ggsave(file.path(OUT, "Fig_validation_CI_coverage.pdf"), p, width = 7.5, height = 5.2)
}

bundle <- list(
  settings = list(replicates = RREP, bootstrap_B = BOOT, bca_replicates = BCA_REP,
                  cores = CORES, equivalence_margin = .30, trajectory_margin = .15,
                  rho = .30, missing_rate = 0),
  scenarios = scenarios,
  failures = failures,
  entity_results = entity_results,
  layer_results = layer_results,
  cell_metrics = cells,
  acceptance = acceptance
)
atomic_save_rds(bundle, file.path(OUT, "confirmatory_validation_bundle.rds"))
capture.output(sessionInfo(), file = file.path(OUT, "sessionInfo.txt"))
message("Confirmatory validation written to: ", normalizePath(OUT, mustWork = FALSE))
