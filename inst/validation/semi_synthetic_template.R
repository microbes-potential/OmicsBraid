# Semi-synthetic validation template
# 1. Load a real matched multi-omics dataset and preprocess each assay appropriately.
# 2. Remove/choose a phenotypically homogeneous background subset.
# 3. Inject known standardized shifts into prespecified pathway scores while preserving
#    the real covariance/noise/missingness background.
# 4. Run OmicsBraid and score recovery of the injected pattern.

library(OmicsBraid)

inject_shift <- function(mat, sample_ids, amount) {
  mat[, sample_ids] <- sweep(mat[, sample_ids, drop = FALSE], 1, amount, "+")
  mat
}

# Recommended injected truths:
# concordant:    +1.0, +1.0, +1.0
# attenuation:   +1.4, +0.8, +0.3
# amplification: +0.3, +0.8, +1.4
# inversion:     +1.0, +0.7, -1.0
# buffering:     +1.0,  0.0,  0.0
# emergence:      0.0,  0.0, +1.0
# null:           0.0,  0.0,  0.0
