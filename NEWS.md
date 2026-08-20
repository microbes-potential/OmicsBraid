# OmicsBraid 0.2.2

## Frozen manuscript release

* Frozen the covariance-aware v0.2.2 statistical core used for manuscript analyses.
* Provides subject-level permutation calibration for the cross-omic omnibus test.
* Provides null-shift matched-subject bootstrap calibration for `Q_omics`.
* Retains analytic confidence intervals as default reporting and BCa layer intervals as sensitivity analysis.
* Retains covariance-aware GLS trajectory inference and hierarchical confirmatory/suggestive braid classification.
* Includes Evidence Forest, Effect Braid, braid heatmap, concordance map, export helpers, and known-truth simulation.
* Adds public-facing GitHub/pkgdown/release documentation without altering the frozen `R/` statistical source.

## Development history

# OmicsBraid 0.2.2

## Robust empirical calibration

- Adds `empirical_omics_tests()` for resampling-calibrated cross-omic inference.
- Global omnibus evidence can be calibrated by matched-subject group-label permutation or by centered matched-subject bootstrap.
- Cross-omic heterogeneity can be calibrated by a raw-data null-shift matched-subject bootstrap (recommended) or an effect-level centered bootstrap under the fitted common-effect null; naive label permutation is deliberately not used for this composite null.
- The permutation omnibus uses a covariance-aware Mahalanobis/Wald statistic calibrated by the joint subject-level permutation distribution. Empirical heterogeneity uses a covariance-aware GLS residual quadratic statistic calibrated by matched-subject bootstrap nulls. Both avoid requiring a chi-square reference distribution under heavy tails.
- Asymptotic p-values are preserved alongside empirical p-values. `run_omics_braid()` reports empirical tests only when explicitly requested; empirical p-values are not made primary unless `empirical_use_as_primary = TRUE`.
- Layer CI workflow now also exposes the already-supported `basic` bootstrap interval.
- Adds a final targeted robust-calibration benchmark comparing asymptotic, permutation, centered-bootstrap, and null-shift-bootstrap Type-I error, inversion power, and analytic/basic/percentile/BCa interval coverage.
- Adds persistent internal-disk checkpointing for the final robust-calibration run (`04_RUN_ROBUST_CALIBRATION.R`).
- The validated v0.1.9/v0.2.1 braid classifier logic is otherwise unchanged.

# OmicsBraid 0.2.1

- I/O-resilience patch for confirmatory validation; statistical algorithms and simulation design are unchanged from v0.2.0.
- Confirmatory checkpoints and high-frequency outputs are now written to persistent internal-disk storage under `~/OmicsBraid_ValidationCache/confirmatory_v020_design`.
- Valid checkpoints from an interrupted v0.2.0 run are imported automatically; incomplete/corrupt RDS files are ignored.
- Checkpoints are validated before reuse and written atomically via temporary-file + rename.
- Final validation outputs are synchronized back to the package `_CONFIRMATORY_VALIDATION_OUTPUT` folder only after the local run completes.

# OmicsBraid 0.2.0

- Froze the v0.1.9 omnibus/GLS/Q/equivalence/trend/classification definitions for confirmatory validation rather than continuing classifier redesign.
- Added `bootstrap_effect_intervals()` with subject-bootstrap percentile, basic, and BCa confidence intervals for layer-specific Hedges' g effects.
- Added `bootstrap_consensus_intervals()` with percentile/basic bootstrap confidence intervals for GLS consensus effects.
- Added end-to-end `ci_method` and `integrated_ci_method` options to `run_omics_braid()` while deliberately retaining analytic SEs and p-values as the inferential basis.
- Preserved analytic intervals alongside bootstrap intervals (`conf_low_analytic`, `conf_high_analytic`) so interval-method sensitivity is auditable.
- Added a targeted confirmatory simulation runner with n/group = 20/40/80/160/320, normal versus heavy-tailed residuals, Monte-Carlo calibration intervals, CI-method comparisons, trend-power curves, equivalence-power curves, covariance-assumption comparators, decisive-classification safety metrics, checkpoint/resume support, and validation figures.
- Added a focused BCa validation subset because BCa requires leave-one-subject-out acceleration and is substantially more computationally expensive.
- Added unit tests ensuring robust intervals are ordered/finite and that changing the displayed CI method does not change analytic p-values or omnibus inference.
- Version 0.2.0 is the confirmatory-validation build motivated by the completed v0.1.9 benchmark, which showed strong core calibration but mild heavy-tail undercoverage for analytic layer CIs.

# OmicsBraid 0.1.9

- Added `test_braid_trend()`, a covariance-aware GLS trajectory test with a prespecified practical slope margin.
- Replaced raw observed-slope attenuation/amplification rules with inferential trajectory states.
- Added hierarchical braid status: confirmed subtype, direction-confirmed broader concordance, no-evidence, unresolved, and insufficient.
- Added `no_detectable_effect` to distinguish failure to reject the joint null from demonstrated practical equivalence (`null_equivalent`).
- Buffering/emergence remain confirmatory only when the required layers pass equivalence testing; a separate `suggestive_pattern` reports effect geometry when precision is insufficient.
- Expanded simulation validation with Monte-Carlo intervals, trend operating characteristics, exact versus hierarchical-family accuracy, null-compatible outcomes, independence-assumption comparators, scenario-failure reporting, and stratification by sample size/correlation/missingness/distribution.

# OmicsBraid 0.1.0

- Initial research implementation.
