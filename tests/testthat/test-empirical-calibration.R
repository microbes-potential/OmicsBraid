test_that("empirical omnibus and heterogeneity calibration return finite p-values", {
  sim <- simulate_braid_data(
    n_per_group = 30,
    patterns = list(joint_null = c(0,0,0), equal_effect = c(.8,.8,.8), inversion = c(.9,.9,-1)),
    rho = .3, seed = 901
  )
  eff <- estimate_effects(sim$data, "group", "Control", "Disease")
  emp <- empirical_omics_tests(
    sim$data, "group", "Control", "Disease",
    effects = eff, B = 99, seed = 902, min_complete = 50,
    omnibus_method = "permutation"
  )
  expect_true(all(c("p_omnibus_empirical", "p_heterogeneity_empirical") %in% names(emp)))
  expect_true(all(is.finite(emp$p_omnibus_empirical)))
  expect_true(all(is.finite(emp$p_heterogeneity_empirical)))
  expect_true(all(emp$p_omnibus_empirical > 0 & emp$p_omnibus_empirical <= 1))
  expect_true(all(emp$p_heterogeneity_empirical > 0 & emp$p_heterogeneity_empirical <= 1))
  expect_equal(emp$omnibus_resamples_used, rep(99L, nrow(emp)))
})

test_that("workflow can report empirical p-values without replacing asymptotic inference", {
  sim <- simulate_braid_data(
    n_per_group = 25,
    patterns = list(concordant_up = c(.9,.9,.9), inversion = c(.9,.9,-1)),
    rho = .25, seed = 911
  )
  fit <- run_omics_braid(
    sim$data, group = "group", reference = "Control", comparison = "Disease",
    bootstrap_B = 100, empirical_tests = TRUE, empirical_B = 99,
    empirical_omnibus_method = "centered_bootstrap",
    empirical_use_as_primary = FALSE,
    pattern_draws = 100, seed = 912
  )
  expect_s3_class(fit, "omics_braid_result")
  expect_true(is.data.frame(fit$empirical_tests))
  expect_true(all(c("p_omnibus_asymptotic", "p_omnibus_empirical",
                    "p_heterogeneity_asymptotic", "p_heterogeneity_empirical") %in%
                  names(fit$integrated)))
  expect_equal(fit$integrated$p_omnibus, fit$integrated$p_omnibus_asymptotic)
  expect_equal(fit$integrated$p_heterogeneity, fit$integrated$p_heterogeneity_asymptotic)
})

test_that("empirical p-values can be selected as primary while preserving asymptotic columns", {
  sim <- simulate_braid_data(
    n_per_group = 25,
    patterns = list(joint_null = c(0,0,0), equal_effect = c(.8,.8,.8)),
    rho = .2, seed = 921
  )
  fit <- run_omics_braid(
    sim$data, group = "group", reference = "Control", comparison = "Disease",
    bootstrap_B = 100, empirical_tests = TRUE, empirical_B = 99,
    empirical_omnibus_method = "centered_bootstrap",
    empirical_use_as_primary = TRUE,
    pattern_draws = 100, seed = 922
  )
  use <- is.finite(fit$integrated$p_omnibus_empirical)
  expect_equal(fit$integrated$p_omnibus[use], fit$integrated$p_omnibus_empirical[use])
  expect_true(all(is.finite(fit$integrated$p_omnibus_asymptotic)))
})
