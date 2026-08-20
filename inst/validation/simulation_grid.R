# OmicsBraid simulation-grid template
# Run after installing OmicsBraid.

library(OmicsBraid)

scenarios <- expand.grid(
  n_per_group = c(20, 40, 80, 150),
  rho = c(0, 0.3, 0.6, 0.85),
  missing_rate = c(0, 0.1, 0.3),
  replicate = seq_len(100),
  KEEP.OUT.ATTRS = FALSE
)

one_run <- function(n_per_group, rho, missing_rate, replicate) {
  sim <- simulate_braid_data(
    n_per_group = n_per_group,
    rho = rho,
    missing_rate = missing_rate,
    seed = 100000 + replicate + n_per_group
  )
  fit <- run_omics_braid(
    sim$data,
    group = "group",
    reference = "Control",
    comparison = "Disease",
    omic_order = c("RNA", "Protein", "Metabolite"),
    bootstrap_B = 300,
    equivalence_margin = 0.30,
    pattern_draws = 1000,
    seed = replicate
  )
  out <- merge(fit$classification, fit$integrated, by = "entity")
  out$n_per_group <- n_per_group
  out$rho <- rho
  out$missing_rate <- missing_rate
  out$replicate <- replicate
  out
}

# For a full benchmark, parallelize this grid and summarize:
# - null integrated-effect type-I error
# - multivariate omnibus type-I error under the all-zero null
# - multivariate omnibus power under concordant and inversion alternatives
# - Q_omics type-I error under truly equal layer effects
# - power of Q_omics under inversion/heterogeneous effects
# - pattern classification accuracy
# - pattern_probability calibration
# - CI coverage
# - sensitivity to independence assumption vs bootstrap covariance
