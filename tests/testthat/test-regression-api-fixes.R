test_that("base-R API regression paths work", {
  effects <- data.frame(
    entity = rep("A", 3),
    omic = c("RNA", "Protein", "Metabolite"),
    effect = c(1, 0.8, -0.9),
    se = rep(0.2, 3)
  )
  x <- integrate_effects(effects)
  expect_equal(nrow(x), 1L)
  expect_true(is.finite(x$Q_omics))

  V <- matrix(c(0.04, 0.01, 0.01, 0.04), 2, 2,
              dimnames = list(c("RNA", "Methylation"), c("RNA", "Methylation")))
  ef <- data.frame(entity = c("A", "A"), omic = c("RNA", "Methylation"),
                   effect = c(1, 0.5), se = c(0.2, 0.2),
                   conf_low = c(0.6, 0.1), conf_high = c(1.4, 0.9))
  z <- orient_omics(ef, list(A = V), c(Methylation = -1))
  expect_equal(z$covariance$A[1, 2], -0.01)

  meta <- data.frame(sample_id = paste0("S", 1:8), group = rep(c("C", "D"), each = 4))
  a <- matrix(rnorm(32), 4, 8, dimnames = list(paste0("G", 1:4), meta$sample_id))
  dat <- omics_braid_data(list(RNA = a), meta)
  mp <- data.frame(omic = "RNA", feature_id = paste0("G", 1:4),
                   pathway = c("P1", "P1", "P1", "P2"))
  sc <- score_pathways(dat, mp, min_features = 2)
  expect_true("P1" %in% rownames(sc$assays$RNA))
})
