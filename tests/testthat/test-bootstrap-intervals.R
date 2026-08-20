test_that("percentile bootstrap intervals are finite and ordered", {
  sim <- simulate_braid_data(
    n_per_group = 20,
    patterns = list(A = c(0.8, 0.8, 0.8)),
    seed = 41
  )
  eff <- estimate_effects(sim$data, "group", "Control", "Disease")
  cv <- bootstrap_effect_covariance(
    sim$data, "group", "Control", "Disease",
    effects = eff, B = 100, seed = 41, min_complete = 50
  )
  ci <- bootstrap_effect_intervals(eff, cv, method = "percentile", min_boot = 50)
  expect_equal(nrow(ci), 3)
  expect_true(all(is.finite(ci$conf_low_boot)))
  expect_true(all(is.finite(ci$conf_high_boot)))
  expect_true(all(ci$conf_low_boot < ci$conf_high_boot))
})

test_that("BCa intervals can be calculated from sample-level data", {
  sim <- simulate_braid_data(
    n_per_group = 15,
    patterns = list(A = c(0.8, 0.8, 0.8)),
    seed = 42
  )
  eff <- estimate_effects(sim$data, "group", "Control", "Disease")
  cv <- bootstrap_effect_covariance(
    sim$data, "group", "Control", "Disease",
    effects = eff, B = 100, seed = 42, min_complete = 50
  )
  ci <- bootstrap_effect_intervals(
    eff, cv, method = "bca", min_boot = 50,
    data = sim$data, group = "group",
    reference = "Control", comparison = "Disease"
  )
  expect_equal(nrow(ci), 3)
  expect_true(all(is.finite(ci$bca_z0)))
  expect_true(all(is.finite(ci$bca_acceleration)))
  expect_true(all(ci$conf_low_boot < ci$conf_high_boot))
})

test_that("workflow can display bootstrap layer and consensus intervals without changing p-values", {
  sim <- simulate_braid_data(
    n_per_group = 20,
    patterns = list(A = c(0.9, 0.8, 0.7)),
    seed = 43
  )
  analytic <- run_omics_braid(
    sim$data, "group", "Control", "Disease",
    bootstrap_B = 100, pattern_draws = 50, seed = 43
  )
  robust <- run_omics_braid(
    sim$data, "group", "Control", "Disease",
    bootstrap_B = 100, ci_method = "percentile",
    integrated_ci_method = "percentile",
    ci_min_boot = 50, pattern_draws = 50, seed = 43
  )
  expect_equal(analytic$effects$p_value, robust$effects$p_value)
  expect_equal(analytic$integrated$p_omnibus, robust$integrated$p_omnibus)
  expect_true(all(c("conf_low_analytic", "conf_low_bootstrap", "ci_method") %in% names(robust$effects)))
  expect_true(all(c("conf_low_analytic", "conf_low_bootstrap", "ci_method") %in% names(robust$integrated)))
})
