test_that("bootstrap covariance is positive semidefinite", {
  sim <- simulate_braid_data(n_per_group = 20, seed = 3)
  ef <- estimate_effects(sim$data, "group", "Control", "Disease")
  cv <- bootstrap_effect_covariance(sim$data, "group", "Control", "Disease",
                                    effects = ef, entities = "concordant_up",
                                    B = 60, min_complete = 20, seed = 3)
  expect_s3_class(cv, "omics_braid_covariance")
  V <- cv$covariance[["concordant_up"]]
  expect_true(all(eigen(V, symmetric = TRUE, only.values = TRUE)$values > -1e-8))
})


test_that("stratified bootstrap records its design", {
  sim <- simulate_braid_data(n_reference = 18, n_comparison = 25, seed = 9)
  ef <- estimate_effects(sim$data, "group", "Control", "Disease")
  cv <- bootstrap_effect_covariance(sim$data, "group", "Control", "Disease", effects = ef,
                                    entities = "attenuation", B = 60, min_complete = 20, seed = 9)
  expect_true(cv$stratified)
})
