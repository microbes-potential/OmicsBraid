# OmicsBraid development checklist

## Before a methods manuscript

- [ ] Replace placeholder package author/maintainer information in DESCRIPTION and CITATION.
- [ ] Run `devtools::document()` with current roxygen2.
- [ ] Run `R CMD build OmicsBraid` and `R CMD check --as-cran` on Linux, macOS, and Windows.
- [ ] Run the full simulation grid with >=1000 replicates for key null scenarios.
- [ ] Compare Hedges' g variance estimators and confidence-interval coverage.
- [ ] Validate generalized Q calibration when covariance is estimated rather than known.
- [ ] Validate GLS trajectory Type-I error/power and sensitivity to trajectory_margin.
- [ ] Validate exact-subtype versus hierarchical-family classification separately.
- [ ] Quantify confirmed versus suggestive buffering/emergence recovery across sample sizes.
- [ ] Quantify the effect of correlation shrinkage and PSD repair.
- [ ] Examine non-normal, heavy-tailed, zero-inflated, and censored distributions.
- [ ] Add cluster/block bootstrap support for repeated measures if longitudinal use is claimed.
- [ ] Add covariate-adjusted estimators only after defining a common estimand across layers.
- [ ] Pre-register/lock braid pattern definitions before real-data discovery benchmarking.
- [ ] Conduct semi-synthetic injection benchmark using at least one real multi-omics cohort.
- [ ] Conduct external real-cohort validation/reproducibility analysis.
- [ ] Compare against qualitative sign overlap, correlation, independent-effect synthesis, and relevant multivariate approaches.
- [ ] Document all sensitivity analyses in a manuscript supplement.

## Scope intentionally excluded from 0.2.0

- raw read/count preprocessing;
- peptide/protein inference;
- metabolite identification or missing-value imputation;
- causal inference;
- longitudinal/repeated-measures effects;
- >2-group omnibus contrasts;
- random-effects distribution over omics layers;
- automatic KEGG/Reactome downloads;
- clinical decision support.
