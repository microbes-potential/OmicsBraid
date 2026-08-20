test_that("GLS trend test distinguishes flat, attenuation, and amplification", {
  effects <- rbind(
    data.frame(entity = "flat", omic = c("RNA","Protein","Metabolite"), effect = c(.8,.8,.8), se = rep(.05,3)),
    data.frame(entity = "att", omic = c("RNA","Protein","Metabolite"), effect = c(1.4,.9,.5), se = rep(.05,3)),
    data.frame(entity = "amp", omic = c("RNA","Protein","Metabolite"), effect = c(.5,.9,1.4), se = rep(.05,3))
  )
  tr <- test_braid_trend(effects, omic_order = c("RNA","Protein","Metabolite"), trajectory_margin = .15)
  expect_equal(tr$trend_state_local[tr$entity == "flat"], "flat")
  expect_equal(tr$trend_state_local[tr$entity == "att"], "attenuation")
  expect_equal(tr$trend_state_local[tr$entity == "amp"], "amplification")
  expect_lt(tr$trend_slope[tr$entity == "att"], 0)
  expect_gt(tr$trend_slope[tr$entity == "amp"], 0)
})

test_that("hierarchical classifier separates null, inversion, buffering and emergence", {
  effects <- rbind(
    data.frame(entity = "inv", omic = c("RNA","Protein","Metabolite"), effect = c(1,.8,-1), se = rep(.05,3)),
    data.frame(entity = "buf", omic = c("RNA","Protein","Metabolite"), effect = c(1.1,0,0), se = rep(.05,3)),
    data.frame(entity = "emg", omic = c("RNA","Protein","Metabolite"), effect = c(0,0,1.1), se = rep(.05,3)),
    data.frame(entity = "none", omic = c("RNA","Protein","Metabolite"), effect = c(.05,-.03,.04), se = rep(.30,3))
  )
  fit <- run_omics_braid_summary(
    effects, omic_order = c("RNA","Protein","Metabolite"),
    equivalence_margin = .30, trajectory_margin = .15, pattern_draws = 100, seed = 1
  )
  cl <- fit$classification
  expect_equal(cl$pattern[cl$entity == "inv"], "inversion")
  expect_equal(cl$pattern[cl$entity == "buf"], "buffering")
  expect_equal(cl$pattern[cl$entity == "emg"], "emergence")
  expect_equal(cl$pattern[cl$entity == "none"], "no_detectable_effect")
})

test_that("same-direction but unresolved trajectory remains broader concordance", {
  effects <- data.frame(
    entity = rep("A",3), omic = c("RNA","Protein","Metabolite"),
    effect = c(1.0,.9,.8), se = rep(.10,3)
  )
  fit <- run_omics_braid_summary(
    effects, omic_order = c("RNA","Protein","Metabolite"),
    equivalence_margin = .30, trajectory_margin = .15, pattern_draws = 100, seed = 2
  )
  expect_equal(fit$classification$pattern, "concordant_increase")
  expect_true(fit$classification$pattern_status %in% c("confirmed", "direction_confirmed"))
})

test_that("mixed-sign trajectories are not assigned attenuation/amplification", {
  effects <- data.frame(
    entity = rep("A",3), omic = c("RNA","Protein","Metabolite"),
    effect = c(1.0,.8,-1.0), se = rep(.08,3)
  )
  tr <- test_braid_trend(effects, omic_order = c("RNA","Protein","Metabolite"))
  expect_equal(tr$trend_state_local, "not_applicable")
})
