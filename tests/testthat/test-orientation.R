test_that("orientation flips effects, intervals, and covariance signs", {
  ef <- data.frame(entity = c("A","A"), omic = c("RNA","Methylation"),
                   effect = c(1,.5), se = c(.2,.2), conf_low = c(.6,.1), conf_high = c(1.4,.9))
  V <- matrix(c(.04,.01,.01,.04),2,2,dimnames=list(c("RNA","Methylation"),c("RNA","Methylation")))
  z <- orient_omics(ef, list(A=V), c(Methylation=-1))
  expect_equal(z$effects$effect[2], -.5)
  expect_equal(z$effects$conf_low[2], -.9)
  expect_equal(z$covariance$A[1,2], -.01)
})
