# OmicsBraid v0.2.0 confirmatory validation metrics

## Calibration

### Omnibus Type-I error
Under the joint-null generating pattern `(0,0,0)`, estimate `Pr(p_omnibus < 0.05)` separately for each sample size and residual distribution. Report Wilson Monte-Carlo 95% intervals. The primary calibration criterion is compatibility of the nominal 0.05 target with the Monte-Carlo interval and absence of systematic drift with n.

### Q_omics Type-I error
Under equal nonzero effects `(0.8,0.8,0.8)`, estimate `Pr(p_heterogeneity < 0.05)` by n/distribution. This tests calibration of the generalized heterogeneity statistic without conflating it with the all-zero null.

### Independence comparator
Repeat omnibus/common-effect/Q calculations with covariance set to independence. This is a deliberately misspecified comparator, not an alternative recommended method.

## Confidence-interval coverage

For every layer effect with known generating effect, record coverage for:

- analytic Hedges-g interval;
- subject-bootstrap percentile interval;
- subject-bootstrap basic interval;
- targeted BCa interval in prespecified heavy-tail subsets.

For equal-effect entities, compare analytic and percentile bootstrap GLS-consensus coverage around the known common effect 0.8.

Bootstrap interval selection must not change analytic SEs or p-values in v0.2.0.

## Ordered-trend operating characteristics

- false meaningful attenuation/amplification calls under equal effects;
- attenuation power by n/distribution;
- amplification power by n/distribution.

## Equivalence operating characteristics

- power to establish practical equivalence for truly zero layers;
- false equivalence rate for true effects beyond the equivalence margin.

The goal is to characterize the sample size needed to confirm buffering/emergence, not to redefine non-significance as equivalence.

## Classification safety

- wrong exact decisive-call rate;
- wrong hierarchical-family decisive-call rate.

Unresolved/suggestive outputs are not counted the same as confidently wrong confirmed outputs.

## Reliability

- scenario failure count/rate;
- finite usable bootstrap-draw counts;
- BCa acceleration/bias-correction diagnostics in the targeted subset.

## Confirmatory simulation grid

Primary grid:

- n/group: 20, 40, 80, 160, 320;
- residuals: normal and heavy-tailed t;
- rho: 0.30;
- imposed missingness: 0;
- Monte-Carlo replicates: 500/cell;
- matched-subject bootstrap draws: 500/simulation;
- BCa subset: first 50 heavy-tail replicates for n=40/80/160.

The wider correlation/missingness sensitivity grid from v0.1.9 is retained as exploratory evidence and is not repeated at the same depth in v0.2.0 unless confirmatory results indicate a need.

# v0.2.2 robust-calibration addendum

## Empirical omnibus calibration

Under the joint-null entity `(0,0,0)`, compare:

- asymptotic covariance-aware chi-square omnibus;
- subject-level permutation omnibus using the joint distribution of `sum(z_k^2)`;
- centered matched-bootstrap omnibus using empirically studentized layer effects.

Type-I error is estimated separately for every sample size and residual distribution with exact binomial Monte-Carlo confidence intervals.

## Empirical heterogeneity calibration

Under equal non-zero effects `(0.8,0.8,0.8)`, compare:

- asymptotic `Q_omics` chi-square p-values;
- centered matched-bootstrap heterogeneity p-values based on effect deviations from the common-effect subspace.

The nominal 0.05 target must be assessed under both normal and heavy-tailed data.

## Robust interval coverage

Compare analytic, percentile, and basic intervals over all simulated layer effects. BCa coverage is targeted to the equal-effect RNA layer in a prespecified replicate subset across every sample size/distribution cell, allowing substantially more BCa validation than v0.2.1 without making the full simulation prohibitively expensive.

## Inversion sensitivity

Because the empirical tests should repair calibration without sacrificing the package's strongest discovery mode, report asymptotic and empirical omnibus/heterogeneity power for the inversion entity `(1.0,0.8,-1.0)`.

## Final acceptance logic

The robust-calibration summary reports PASS/REVIEW for:

1. permutation omnibus Type-I calibration;
2. centered-bootstrap heterogeneity Type-I calibration;
3. heavy-tail BCa interval coverage;
4. empirical inversion sensitivity at n >= 40;
5. zero scenario-execution failures.

These criteria are validation gates, not claims that every finite-simulation point estimate must equal its nominal target exactly.

The v0.2.2 robust benchmark reports both the effect-level centered-bootstrap heterogeneity calibration and the raw-data **null-shift bootstrap**. The final acceptance gate is based on the null-shift method because it explicitly imposes equal standardized effects while preserving residual tail behavior.
