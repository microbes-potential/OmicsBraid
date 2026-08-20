test_that("simulation constructs valid multi-omics data", {
  sim <- simulate_braid_data(n_per_group = 20, seed = 1)
  expect_s3_class(sim$data, "omics_braid_data")
  expect_equal(length(sim$data$assays), 3)
  expect_true(all(c("Control", "Disease") %in% sim$data$metadata$group))
})

test_that("effects recover intended directions", {
  sim <- simulate_braid_data(n_per_group = 120, rho = 0.2, seed = 2)
  ef <- estimate_effects(sim$data, "group", "Control", "Disease")
  inv <- ef[ef$entity == "inversion", ]
  expect_true(inv$effect[inv$omic == "RNA"] > 0)
  expect_true(inv$effect[inv$omic == "Metabolite"] < 0)
})
