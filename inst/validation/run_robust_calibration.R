# OmicsBraid v0.2.2 final robust-calibration validation
#
# Targets only the remaining issues isolated by the v0.2.1 confirmatory study:
#   1) chi-square omnibus calibration under heavy tails;
#   2) chi-square Q_omics calibration under heavy tails;
#   3) empirical permutation/centered-bootstrap calibration;
#   4) analytic/basic/percentile/BCa layer-CI coverage;
#   5) inversion power after empirical calibration.
#
# Example:
# Rscript run_robust_calibration.R --replicates=500 --B-cov=300 --B-emp=199 --bca-replicates=100 --cores=2 --chunk=20 --out=robust_results

validation_lib <- Sys.getenv("OMICSBRAID_VALIDATION_LIB", unset = "")
if (nzchar(validation_lib)) {
  if (!dir.exists(validation_lib)) stop("OMICSBRAID_VALIDATION_LIB does not exist: ", validation_lib, call. = FALSE)
  .libPaths(unique(c(normalizePath(validation_lib, mustWork = TRUE), .libPaths())))
}
suppressPackageStartupMessages(library(OmicsBraid))
expected_version <- Sys.getenv("OMICSBRAID_EXPECTED_VERSION", unset = "")
if (nzchar(expected_version)) {
  iv <- as.character(utils::packageVersion("OmicsBraid"))
  if (!identical(iv, expected_version)) stop("Loaded OmicsBraid ", iv, " but expected ", expected_version, call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", key, "="), "", hit[1])
}
RREP <- as.integer(get_arg("replicates", "500"))
B_COV <- as.integer(get_arg("B-cov", "300"))
B_EMP <- as.integer(get_arg("B-emp", "199"))
BCA_REP <- as.integer(get_arg("bca-replicates", "100"))
CORES <- as.integer(get_arg("cores", "2"))
CHUNK <- as.integer(get_arg("chunk", "20"))
OUT <- get_arg("out", "OmicsBraid_robust_calibration")
if (!is.finite(RREP) || RREP < 100L) stop("--replicates must be >= 100", call. = FALSE)
if (!is.finite(B_COV) || B_COV < 100L) stop("--B-cov must be >= 100", call. = FALSE)
if (!is.finite(B_EMP) || B_EMP < 99L) stop("--B-emp must be >= 99", call. = FALSE)
if (!is.finite(BCA_REP) || BCA_REP < 0L) BCA_REP <- 0L
if (!is.finite(CORES) || CORES < 1L) CORES <- 1L
if (!is.finite(CHUNK) || CHUNK < 1L) CHUNK <- 20L
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
checkpoint_dir <- file.path(OUT, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

checkpoint_is_valid <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0) return(FALSE)
  z <- try(readRDS(path), silent = TRUE)
  !inherits(z, "try-error") && is.list(z) && length(z) > 0L
}
atomic_save_rds <- function(object, path, attempts = 5L) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  for (attempt in seq_len(attempts)) {
    tmp <- tempfile(pattern = paste0(".", basename(path), ".tmp_"), tmpdir = dirname(path), fileext = ".rds")
    ok <- tryCatch({
      saveRDS(object, tmp, compress = FALSE)
      if (!file.exists(tmp) || file.info(tmp)$size <= 0) stop("temporary checkpoint missing")
      if (file.exists(path)) unlink(path, force = TRUE)
      if (!file.rename(tmp, path)) stop("checkpoint rename failed")
      checkpoint_is_valid(path)
    }, error = function(e) FALSE)
    unlink(tmp, force = TRUE)
    if (isTRUE(ok)) return(invisible(TRUE))
    Sys.sleep(.25 * attempt)
  }
  stop("Unable to write checkpoint: ", path, call. = FALSE)
}

omics <- c("RNA", "Protein", "Metabolite")
truth_patterns <- list(
  joint_null = c(0, 0, 0),
  equal_effect = c(.8, .8, .8),
  inversion = c(1.0, .8, -1.0)
)
scenarios <- expand.grid(
  n_per_group = c(20L, 40L, 80L, 160L, 320L),
  residual_distribution = c("normal", "t"),
  replicate = seq_len(RREP),
  stringsAsFactors = FALSE,
  KEEP.OUT.ATTRS = FALSE
)
scenarios$rho <- .30
scenarios$scenario_id <- seq_len(nrow(scenarios))

seed_for <- function(z) {
  52000000L + as.integer(z$replicate) +
    10000L * as.integer(z$n_per_group) +
    1000000L * match(z$residual_distribution, c("normal", "t"))
}

one_run <- function(z) {
  seed <- seed_for(z)
  sim <- simulate_braid_data(
    n_per_group = z$n_per_group,
    omics = omics, patterns = truth_patterns, rho = z$rho,
    missing_rate = 0, residual_distribution = z$residual_distribution,
    seed = seed
  )
  eff <- estimate_effects(sim$data, "group", "Control", "Disease", conf_level = .95)
  cv <- bootstrap_effect_covariance(
    sim$data, "group", "Control", "Disease", effects = eff,
    B = B_COV, seed = seed + 1001L, shrinkage = .05,
    min_complete = max(50L, floor(B_COV * .50)), stratified = TRUE
  )
  asym <- integrate_effects(eff, covariance = cv, conf_level = .95)
  emp <- empirical_omics_tests(
    sim$data, "group", "Control", "Disease", effects = eff,
    B = B_EMP, seed = seed + 2003L,
    min_complete = max(50L, floor(B_EMP * .50)),
    omnibus_method = "permutation",
    heterogeneity_method = "null_shift_bootstrap"
  )

  pct <- bootstrap_effect_intervals(eff, cv, method = "percentile", conf_level = .95,
                                    min_boot = max(50L, floor(B_COV * .50)))
  basic <- bootstrap_effect_intervals(eff, cv, method = "basic", conf_level = .95,
                                      min_boot = max(50L, floor(B_COV * .50)))

  # BCa is targeted to the equal-effect RNA layer only. That is sufficient for
  # direct coverage calibration while avoiding unnecessary jackknife work.
  bca <- NULL
  if (BCA_REP > 0L && z$replicate <= BCA_REP) {
    one_eff <- eff[eff$entity == "equal_effect" & eff$omic == "RNA", , drop = FALSE]
    bca <- bootstrap_effect_intervals(
      one_eff, cv, method = "bca", conf_level = .95,
      min_boot = max(50L, floor(B_COV * .50)), data = sim$data,
      group = "group", reference = "Control", comparison = "Disease"
    )
  }

  entity <- merge(asym, emp, by = "entity", all.x = TRUE, sort = FALSE)
  entity$n_per_group <- z$n_per_group
  entity$residual_distribution <- z$residual_distribution
  entity$replicate <- z$replicate
  entity$scenario_id <- z$scenario_id

  layer <- merge(eff, sim$truth, by = c("entity", "omic"), all.x = TRUE, sort = FALSE)
  key <- paste(layer$entity, layer$omic, sep = "\r")
  kp <- paste(pct$entity, pct$omic, sep = "\r")
  kb <- paste(basic$entity, basic$omic, sep = "\r")
  ip <- match(key, kp); ib <- match(key, kb)
  layer$analytic_covered95 <- layer$conf_low <= layer$true_effect & layer$conf_high >= layer$true_effect
  layer$percentile_low <- pct$conf_low_boot[ip]
  layer$percentile_high <- pct$conf_high_boot[ip]
  layer$percentile_covered95 <- layer$percentile_low <= layer$true_effect & layer$percentile_high >= layer$true_effect
  layer$basic_low <- basic$conf_low_boot[ib]
  layer$basic_high <- basic$conf_high_boot[ib]
  layer$basic_covered95 <- layer$basic_low <= layer$true_effect & layer$basic_high >= layer$true_effect
  layer$bca_low <- NA_real_; layer$bca_high <- NA_real_; layer$bca_covered95 <- NA
  if (!is.null(bca) && nrow(bca)) {
    hit <- layer$entity == "equal_effect" & layer$omic == "RNA"
    layer$bca_low[hit] <- bca$conf_low_boot[1]
    layer$bca_high[hit] <- bca$conf_high_boot[1]
    layer$bca_covered95[hit] <- bca$conf_low_boot[1] <= .8 && bca$conf_high_boot[1] >= .8
  }
  layer$n_per_group <- z$n_per_group
  layer$residual_distribution <- z$residual_distribution
  layer$replicate <- z$replicate
  layer$scenario_id <- z$scenario_id
  list(entity = entity, layer = layer)
}

run_safe <- function(i) {
  z <- scenarios[i, , drop = FALSE]
  tryCatch(
    list(ok = TRUE, result = one_run(z), scenario_id = i, error = NA_character_),
    error = function(e) list(ok = FALSE, result = NULL, scenario_id = i,
                             error = conditionMessage(e))
  )
}

ids <- scenarios$scenario_id
chunks <- split(ids, ceiling(seq_along(ids) / CHUNK))
for (cc in seq_along(chunks)) {
  f <- file.path(checkpoint_dir, sprintf("chunk_%04d.rds", cc))
  if (checkpoint_is_valid(f)) next
  if (file.exists(f)) unlink(f, force = TRUE)
  idx <- chunks[[cc]]
  if (.Platform$OS.type == "unix" && CORES > 1L) {
    zz <- parallel::mclapply(idx, run_safe, mc.cores = CORES, mc.preschedule = FALSE)
  } else {
    zz <- lapply(idx, run_safe)
  }
  atomic_save_rds(zz, f)
  message("Completed robust-calibration chunk ", cc, " / ", length(chunks),
          " (scenario rows through ", max(idx), " / ", nrow(scenarios), ")")
}

ck <- sort(list.files(checkpoint_dir, pattern = "^chunk_[0-9]+\\.rds$", full.names = TRUE))
ck <- ck[vapply(ck, checkpoint_is_valid, logical(1))]
allruns <- unlist(lapply(ck, readRDS), recursive = FALSE)
failures <- do.call(rbind, lapply(allruns[!vapply(allruns, `[[`, logical(1), "ok")], function(z) {
  data.frame(scenario_id = z$scenario_id, error = z$error, stringsAsFactors = FALSE)
}))
if (is.null(failures)) failures <- data.frame(scenario_id = integer(), error = character(), stringsAsFactors = FALSE)
ok <- allruns[vapply(allruns, `[[`, logical(1), "ok")]
entity_results <- do.call(rbind, lapply(ok, function(z) z$result$entity))
layer_results <- do.call(rbind, lapply(ok, function(z) z$result$layer))
rownames(entity_results) <- NULL; rownames(layer_results) <- NULL
utils::write.csv(failures, file.path(OUT, "scenario_failures.csv"), row.names = FALSE)
utils::write.csv(entity_results, file.path(OUT, "entity_level_robust_calibration.csv"), row.names = FALSE)
utils::write.csv(layer_results, file.path(OUT, "layer_level_robust_calibration.csv"), row.names = FALSE)

binom_ci <- function(x, n, conf = .95) {
  if (!n) return(c(NA_real_, NA_real_))
  a <- 1 - conf
  lo <- if (x == 0) 0 else stats::qbeta(a/2, x, n - x + 1)
  hi <- if (x == n) 1 else stats::qbeta(1 - a/2, x + 1, n - x)
  c(lo, hi)
}
metric <- function(name, v, target = NA_real_) {
  v <- v[!is.na(v)]
  n <- length(v); x <- sum(v)
  val <- if (n) mean(v) else NA_real_
  ci <- if (n) binom_ci(x, n) else c(NA_real_, NA_real_)
  data.frame(metric = name, value = val, n_evaluable = n,
             mc_ci_low = ci[1], mc_ci_high = ci[2], target = target,
             target_inside_mc_ci = if (is.finite(target) && n) ci[1] <= target && ci[2] >= target else NA,
             stringsAsFactors = FALSE)
}

cell_summary <- function(n, dist) {
  e <- entity_results[entity_results$n_per_group == n & entity_results$residual_distribution == dist, , drop = FALSE]
  l <- layer_results[layer_results$n_per_group == n & layer_results$residual_distribution == dist, , drop = FALSE]
  null <- e$entity == "joint_null"
  equal <- e$entity == "equal_effect"
  inv <- e$entity == "inversion"
  bca_rows <- l$entity == "equal_effect" & l$omic == "RNA" & !is.na(l$bca_covered95)
  rr <- rbind(
    metric("omnibus_asymptotic_type_I", e$p_omnibus[null] < .05, .05),
    metric("omnibus_permutation_type_I", e$p_omnibus_permutation[null] <= .05, .05),
    metric("omnibus_bootstrap_type_I", e$p_omnibus_centered_bootstrap[null] <= .05, .05),
    metric("Q_asymptotic_type_I", e$p_heterogeneity[equal] < .05, .05),
    metric("Q_centered_bootstrap_type_I", e$p_heterogeneity_centered_bootstrap[equal] <= .05, .05),
    metric("Q_null_shift_bootstrap_type_I", e$p_heterogeneity_null_shift_bootstrap[equal] <= .05, .05),
    metric("omnibus_asymptotic_power_inversion", e$p_omnibus[inv] < .05),
    metric("omnibus_permutation_power_inversion", e$p_omnibus_permutation[inv] <= .05),
    metric("omnibus_bootstrap_power_inversion", e$p_omnibus_centered_bootstrap[inv] <= .05),
    metric("Q_asymptotic_power_inversion", e$p_heterogeneity[inv] < .05),
    metric("Q_centered_bootstrap_power_inversion", e$p_heterogeneity_centered_bootstrap[inv] <= .05),
    metric("Q_null_shift_bootstrap_power_inversion", e$p_heterogeneity_null_shift_bootstrap[inv] <= .05),
    metric("analytic_layer_CI_coverage", l$analytic_covered95, .95),
    metric("percentile_layer_CI_coverage", l$percentile_covered95, .95),
    metric("basic_layer_CI_coverage", l$basic_covered95, .95),
    metric("BCa_equal_RNA_CI_coverage", l$bca_covered95[bca_rows], .95)
  )
  rr$n_per_group <- n; rr$residual_distribution <- dist
  rr[, c("n_per_group","residual_distribution","metric","value","n_evaluable","mc_ci_low","mc_ci_high","target","target_inside_mc_ci")]
}
metrics <- do.call(rbind, lapply(c(20L,40L,80L,160L,320L), function(n) {
  do.call(rbind, lapply(c("normal","t"), function(d) cell_summary(n,d)))
}))
utils::write.csv(metrics, file.path(OUT, "robust_metrics_by_sample_size_distribution.csv"), row.names = FALSE)

pick <- function(x) metrics[metrics$metric %in% x, , drop = FALSE]
utils::write.csv(pick(c("omnibus_asymptotic_type_I","omnibus_permutation_type_I","omnibus_bootstrap_type_I",
                        "Q_asymptotic_type_I","Q_centered_bootstrap_type_I","Q_null_shift_bootstrap_type_I")),
                 file.path(OUT, "robust_type1_calibration.csv"), row.names = FALSE)
utils::write.csv(pick(c("omnibus_asymptotic_power_inversion","omnibus_permutation_power_inversion",
                        "omnibus_bootstrap_power_inversion","Q_asymptotic_power_inversion",
                        "Q_centered_bootstrap_power_inversion","Q_null_shift_bootstrap_power_inversion")),
                 file.path(OUT, "robust_inversion_power.csv"), row.names = FALSE)
utils::write.csv(pick(c("analytic_layer_CI_coverage","percentile_layer_CI_coverage",
                        "basic_layer_CI_coverage","BCa_equal_RNA_CI_coverage")),
                 file.path(OUT, "robust_CI_coverage.csv"), row.names = FALSE)

# Aggregate CI coverage by distribution for easy comparison.
ci_dist <- do.call(rbind, lapply(c("normal","t"), function(d) {
  l <- layer_results[layer_results$residual_distribution == d, , drop = FALSE]
  bca_rows <- l$entity == "equal_effect" & l$omic == "RNA" & !is.na(l$bca_covered95)
  rr <- rbind(
    metric("analytic_layer_CI_coverage", l$analytic_covered95, .95),
    metric("percentile_layer_CI_coverage", l$percentile_covered95, .95),
    metric("basic_layer_CI_coverage", l$basic_covered95, .95),
    metric("BCa_equal_RNA_CI_coverage", l$bca_covered95[bca_rows], .95)
  )
  rr$residual_distribution <- d
  rr[, c("residual_distribution","metric","value","n_evaluable","mc_ci_low","mc_ci_high","target","target_inside_mc_ci")]
}))
utils::write.csv(ci_dist, file.path(OUT, "robust_CI_coverage_by_distribution.csv"), row.names = FALSE)

# Acceptance criteria: the robust methods must repair the specific heavy-tail
# weaknesses without sacrificing inversion sensitivity.
emp_omni <- metrics[metrics$metric == "omnibus_permutation_type_I", , drop = FALSE]
emp_q <- metrics[metrics$metric == "Q_null_shift_bootstrap_type_I", , drop = FALSE]
bca_t <- metrics[metrics$metric == "BCa_equal_RNA_CI_coverage" & metrics$residual_distribution == "t" & metrics$n_evaluable > 0, , drop = FALSE]
pow_inv <- metrics[metrics$metric %in% c("omnibus_permutation_power_inversion","Q_null_shift_bootstrap_power_inversion") & metrics$n_per_group >= 40, , drop = FALSE]
acceptance <- rbind(
  data.frame(criterion = "Permutation omnibus Type-I calibration",
             result = if (nrow(emp_omni) && all(emp_omni$target_inside_mc_ci %in% TRUE)) "PASS" else "REVIEW",
             rule = "Nominal 0.05 lies inside every n x distribution Monte Carlo 95% interval", stringsAsFactors = FALSE),
  data.frame(criterion = "Null-shift-bootstrap Q Type-I calibration",
             result = if (nrow(emp_q) && all(emp_q$target_inside_mc_ci %in% TRUE)) "PASS" else "REVIEW",
             rule = "Nominal 0.05 lies inside every n x distribution Monte Carlo 95% interval", stringsAsFactors = FALSE),
  data.frame(criterion = "Heavy-tail BCa layer-CI coverage",
             result = if (nrow(bca_t) && all(bca_t$target_inside_mc_ci %in% TRUE)) "PASS" else "REVIEW",
             rule = "Nominal 0.95 lies inside each evaluated heavy-tail sample-size Monte Carlo interval", stringsAsFactors = FALSE),
  data.frame(criterion = "Empirical inversion sensitivity",
             result = if (nrow(pow_inv) && min(pow_inv$value, na.rm = TRUE) >= .90) "PASS" else "REVIEW",
             rule = "Permutation omnibus and bootstrap Q power >= 0.90 for n >= 40", stringsAsFactors = FALSE),
  data.frame(criterion = "Scenario execution reliability",
             result = if (nrow(failures) == 0L) "PASS" else "REVIEW",
             rule = "No failed scenario rows", stringsAsFactors = FALSE)
)
utils::write.csv(acceptance, file.path(OUT, "robust_acceptance_summary.csv"), row.names = FALSE)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_metric <- function(metric_names, ylab, target, filename) {
    d <- metrics[metrics$metric %in% metric_names & is.finite(metrics$value), , drop = FALSE]
    p <- ggplot2::ggplot(d, ggplot2::aes(x = n_per_group, y = value,
                                         linetype = metric, shape = residual_distribution,
                                         group = interaction(metric, residual_distribution))) +
      ggplot2::geom_line() + ggplot2::geom_point() +
      ggplot2::geom_hline(yintercept = target, linetype = 3) +
      ggplot2::scale_x_continuous(breaks = c(20,40,80,160,320)) +
      ggplot2::labs(x = "Sample size per group", y = ylab,
                    linetype = "Method", shape = "Residual distribution") +
      ggplot2::theme_classic()
    ggplot2::ggsave(file.path(OUT, filename), p, width = 8, height = 5.2)
  }
  plot_metric(c("omnibus_asymptotic_type_I","omnibus_permutation_type_I","omnibus_bootstrap_type_I"),
              "Omnibus Type-I error", .05, "Fig_robust_omnibus_type1.pdf")
  plot_metric(c("Q_asymptotic_type_I","Q_centered_bootstrap_type_I","Q_null_shift_bootstrap_type_I"),
              "Qomics Type-I error", .05, "Fig_robust_Q_type1.pdf")
  plot_metric(c("analytic_layer_CI_coverage","percentile_layer_CI_coverage","basic_layer_CI_coverage","BCa_equal_RNA_CI_coverage"),
              "95% CI coverage", .95, "Fig_robust_CI_coverage.pdf")
}

bundle <- list(settings = list(replicates = RREP, B_cov = B_COV, B_emp = B_EMP,
                               bca_replicates = BCA_REP, rho = .30,
                               sample_sizes = c(20,40,80,160,320),
                               distributions = c("normal","t")),
               scenarios = scenarios, failures = failures,
               entity_results = entity_results, layer_results = layer_results,
               metrics = metrics, acceptance = acceptance)
atomic_save_rds(bundle, file.path(OUT, "robust_calibration_bundle.rds"))
capture.output(sessionInfo(), file = file.path(OUT, "sessionInfo.txt"))
message("Robust calibration written to: ", normalizePath(OUT, mustWork = FALSE))
