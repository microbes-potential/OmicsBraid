test_that("summary workflow labels covariance mode", {
  effects <- data.frame(entity = rep("A",3), omic = c("RNA","Protein","Metabolite"),
                        effect = c(.8,.7,.6), se = rep(.2,3))
  fit <- run_omics_braid_summary(effects, omic_order = c("RNA","Protein","Metabolite"), pattern_draws = 100)
  expect_equal(fit$integrated$covariance_mode, "independence")
})

test_that("partial supplied covariance retains all effects", {
  effects <- data.frame(entity = rep("A",3), omic = c("RNA","Protein","Metabolite"),
                        effect = c(.8,.7,.6), se = rep(.2,3))
  V <- matrix(c(.04,.01,.01,.04), 2, 2, dimnames = list(c("RNA","Protein"), c("RNA","Protein")))
  x <- integrate_effects(effects, covariance = list(A = V))
  expect_equal(x$n_omics, 3)
  expect_equal(x$covariance_mode, "partial_supplied")
})
