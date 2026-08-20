test_that("integration returns heterogeneity diagnostics", {
  effects <- data.frame(entity = rep("A", 3), omic = c("RNA","Protein","Metabolite"),
                        effect = c(1, .9, .8), se = c(.2,.2,.2))
  x <- integrate_effects(effects)
  expect_equal(nrow(x), 1)
  expect_true(all(c("W_omnibus", "p_omnibus", "Q_omics", "I2_omics", "integrated_effect") %in% names(x)))
  expect_true(x$integrated_effect > 0)
})

test_that("equivalence distinguishes negligible effects", {
  effects <- data.frame(entity = c("A","B"), omic = c("RNA","RNA"),
                        effect = c(.01, 1.0), se = c(.05,.1),
                        p_value = c(.84, 1e-10), p_adj = c(.84, 1e-10), df = c(100,100))
  eq <- test_equivalence(effects, margin = .3, p_adjust = "none")
  expect_true(eq$equivalent[1])
  expect_false(eq$equivalent[2])
})

test_that("deterministic braid labels inversion", {
  effects <- data.frame(entity = rep("A", 3), omic = c("RNA","Protein","Metabolite"),
                        effect = c(1,.8,-1), se = c(.1,.1,.1),
                        p_value = c(1e-10,1e-8,1e-10), p_adj = c(1e-10,1e-8,1e-10), df = rep(100,3))
  eq <- test_equivalence(effects, margin = .3, p_adjust = "none")
  cl <- classify_braids(eq, c("RNA","Protein","Metabolite"))
  expect_equal(cl$pattern, "inversion")
})


test_that("omnibus detects strong inversion even when consensus cancels", {
  effects <- data.frame(entity = rep("A", 3), omic = c("RNA","Protein","Metabolite"),
                        effect = c(1, 1, -1), se = c(.1,.1,.1))
  x <- integrate_effects(effects)
  expect_true(abs(x$integrated_effect) < 0.5)
  expect_true(x$p_omnibus < 1e-6)
  expect_true(x$p_heterogeneity < 1e-6)
})

test_that("local braid states are separated from multiplicity-adjusted states", {
  effects <- data.frame(
    entity = c("A", "A", paste0("N", 1:20)),
    omic = c("RNA", "Protein", rep("RNA", 20)),
    effect = c(1.0, 0.01, rep(0, 20)),
    se = c(0.10, 0.05, rep(0.5, 20)),
    p_value = c(1e-12, 0.84, rep(1, 20))
  )
  # Avoid duplicate entity x omic rows while retaining many unrelated tests.
  effects$entity[3:22] <- paste0("N", seq_len(20))
  eq <- test_equivalence(effects, margin = 0.3, state_basis = "local")
  a <- eq[eq$entity == "A", ]
  expect_equal(a$state[a$omic == "RNA"], "positive")
  expect_equal(a$state[a$omic == "Protein"], "equivalent")
  expect_true(all(c("state_local", "state_adjusted", "p_difference", "p_difference_adj") %in% names(eq)))
})

test_that("classification reports uncertainty support without replacing deterministic label", {
  effects <- data.frame(entity = rep("A", 3), omic = c("RNA", "Protein", "Metabolite"),
                        effect = c(1, .8, -1), se = c(.1,.1,.1))
  fit <- run_omics_braid_summary(effects, omic_order = c("RNA","Protein","Metabolite"),
                                 equivalence_margin = .3, pattern_draws = 200, seed = 3)
  expect_equal(fit$classification$pattern, "inversion")
  expect_true(all(c("uncertainty_mode_pattern", "uncertainty_mode_support",
                    "deterministic_pattern_support", "pattern_entropy",
                    "classification_stability") %in% names(fit$classification)))
})

test_that("directional alignment is evidence-qualified under the joint null", {
  effects <- data.frame(entity = rep(c("Signal", "Null"), each = 3),
                        omic = rep(c("RNA","Protein","Metabolite"), 2),
                        effect = c(1.0,.9,.8, .05,.04,.03),
                        se = rep(.15, 6))
  x <- integrate_effects(effects)
  sig <- x[x$entity == "Signal", ]
  nul <- x[x$entity == "Null", ]
  expect_true(sig$has_omnibus_evidence)
  expect_true(is.finite(sig$evidence_qualified_direction_agreement))
  expect_false(nul$has_omnibus_evidence)
  expect_true(is.na(nul$evidence_qualified_direction_agreement))
})
