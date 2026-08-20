test_that("plots are ggplot objects", {
  effects <- data.frame(entity = rep("A",3), omic = c("RNA","Protein","Metabolite"),
                        effect = c(1,.8,-.7), se = c(.2,.2,.2),
                        conf_low = c(.6,.4,-1.1), conf_high = c(1.4,1.2,-.3),
                        p_value = c(1e-5,1e-4,1e-3), p_adj = c(1e-5,1e-4,1e-3), df = rep(50,3))
  fit <- run_omics_braid_summary(effects, omic_order = c("RNA","Protein","Metabolite"), pattern_draws = 100)
  expect_s3_class(plot_evidence_forest(fit, "A"), "ggplot")
  expect_s3_class(plot_effect_braid(fit, "A"), "ggplot")
})
