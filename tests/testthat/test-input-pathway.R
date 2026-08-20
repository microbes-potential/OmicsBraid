test_that("partially overlapping assays are accepted", {
  meta <- data.frame(sample_id = paste0("S", 1:6), group = rep(c("C","D"), each = 3))
  a <- matrix(rnorm(12), 2, 6, dimnames = list(c("A","B"), meta$sample_id))
  b <- matrix(rnorm(8), 2, 4, dimnames = list(c("A","B"), meta$sample_id[1:4]))
  x <- omics_braid_data(list(RNA = a, Protein = b), meta)
  expect_s3_class(x, "omics_braid_data")
  expect_equal(ncol(x$assays$Protein), 4)
})

test_that("pathway scoring returns pathway matrices", {
  meta <- data.frame(sample_id = paste0("S", 1:8), group = rep(c("C","D"), each = 4))
  a <- matrix(rnorm(32), 4, 8, dimnames = list(paste0("G",1:4), meta$sample_id))
  x <- omics_braid_data(list(RNA = a), meta)
  mp <- data.frame(omic = "RNA", feature_id = paste0("G",1:4),
                   pathway = c("P1","P1","P1","P2"))
  y <- score_pathways(x, mp, min_features = 2)
  expect_true("P1" %in% rownames(y$assays$RNA))
  expect_false("P2" %in% rownames(y$assays$RNA))
})


test_that("duplicate entity-omic rows are rejected", {
  effects <- data.frame(entity = c("A","A"), omic = c("RNA","RNA"), effect = c(.2,.3), se = c(.1,.1))
  expect_error(integrate_effects(effects), "must occur once")
})

test_that("infinite assay values are rejected", {
  meta <- data.frame(sample_id = c("S1","S2","S3","S4"), group = c("C","C","D","D"))
  a <- matrix(c(1,2,3,Inf), 1, 4, dimnames = list("A", meta$sample_id))
  expect_error(omics_braid_data(list(RNA = a), meta), "infinite values")
})

test_that("entity harmonization maps IDs and rejects ambiguous collapse by default", {
  meta <- data.frame(sample_id = paste0("S",1:6), group = rep(c("C","D"), each = 3))
  rna <- matrix(rnorm(12), 2, 6, dimnames = list(c("ENSG1","ENSG2"), meta$sample_id))
  pro <- matrix(rnorm(12), 2, 6, dimnames = list(c("P1","P2"), meta$sample_id))
  x <- omics_braid_data(list(RNA=rna, Protein=pro), meta)
  mp <- data.frame(omic=c("RNA","RNA","Protein","Protein"),
                   feature_id=c("ENSG1","ENSG2","P1","P2"),
                   entity=c("G1","G2","G1","G2"))
  y <- harmonize_entities(x, mp)
  expect_equal(rownames(y$assays$RNA), c("G1","G2"))
  expect_equal(rownames(y$assays$Protein), c("G1","G2"))
})
