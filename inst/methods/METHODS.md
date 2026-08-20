# OmicsBraid statistical specification (v0.2.0)

## 1. Target estimand

For biological entity/pathway `j` and ordered omic layer `k`, OmicsBraid estimates a standardized contrast

`delta_jk = E(X_jk | comparison) - E(X_jk | reference)`

on a within-omic pooled-standard-deviation scale. Version 0.2.0 implements the small-sample corrected standardized mean difference (Hedges' g).

The primary object of inference is the **vector** `delta_j = (delta_j1, ..., delta_jK)`, not a latent multi-omics factor.

## 2. Layer-specific effect

Let `d` denote Cohen's standardized mean difference using the pooled SD. OmicsBraid applies the gamma-function small-sample correction `J(df)` to obtain `g = J(df)d`.

The analytic variance currently uses the common large-sample approximation

`Var(g) ~= (n0+n1)/(n0*n1) + g^2/(2*(n0+n1))`.

Alternative SMD variance estimators remain a sensitivity-analysis target for a manuscript-grade release.

## 3. Cross-omic covariance

For matched-subject studies, layer-specific effects are statistically dependent. OmicsBraid resamples **subjects**, preserving all available modalities from each sampled subject. By default the bootstrap is stratified by comparison group, preserving the observed group sizes while resampling biological units within group.

For each entity, the bootstrap estimates an effect correlation matrix `R_boot`. The working covariance is

`V = D_se R_boot D_se`,

where `D_se` contains analytic marginal standard errors. `R_boot` is shrunk toward the identity and projected to positive-semidefinite form when necessary. Diagnostics are retained for audit.

## 4. Effect orientation

Concordance has meaning only after effect directions share a prespecified interpretation. OmicsBraid supports an optional omic-level multiplier (+1/-1). If `D` is the diagonal orientation matrix,

`y* = D y` and `V* = D V D`.

Orientation must be specified from biological semantics before examining concordance and must never be optimized to maximize agreement.

## 5. Multivariate omnibus evidence

Before pooling, OmicsBraid tests

`H0: delta_1 = ... = delta_K = 0`

with

`W_omnibus = y' V^-1 y`,

using a large-sample chi-square reference with `K` degrees of freedom.

The omnibus test is intentionally distinct from the consensus effect. A strong inversion such as `(+, +, -)` can have a consensus near zero but overwhelmingly reject the all-zero null.

## 6. GLS consensus effect

For effect vector `y` and covariance `V`,

`theta_hat = (1' V^-1 y) / (1' V^-1 1)`

and

`SE(theta_hat) = sqrt[1 / (1' V^-1 1)]`.

This is a descriptive/common-effect summary, not a causal parameter. When heterogeneity is strong, the layer effects and braid structure should be emphasized over the pooled estimate.

## 7. Cross-omic heterogeneity

OmicsBraid computes

`Q_omics = (y - theta_hat 1)' V^-1 (y - theta_hat 1)`.

Under a correctly specified common-effect null and known covariance, the large-sample reference is chi-square with `K-1` degrees of freedom.

The descriptive index

`I2_omics = max(0, (Q_omics - (K-1))/Q_omics) * 100%`

must not be described as literal meta-analytic between-study heterogeneity. Omics layers are not studies. No universal low/moderate/high cutoffs are imposed in v0.2.0.

## 8. Robust confidence-interval layer

The analytic Hedges-g standard error and its normal-approximation confidence interval remain the inferential reference used by the omnibus, equivalence, and trajectory procedures in v0.2.0. Confirmatory validation of v0.1.9 showed mild analytic interval undercoverage under heavy-tailed residuals, motivating an additional uncertainty-reporting layer rather than an unvalidated replacement of the tests.

Given matched subject-bootstrap effect draws `g^(b)`, OmicsBraid can report:

- **percentile interval:** empirical `alpha/2` and `1-alpha/2` quantiles of `g^(b)`;
- **basic interval:** reflection of bootstrap quantiles around the observed effect;
- **BCa interval:** bias-corrected and accelerated quantiles, with bias correction estimated from the bootstrap distribution and acceleration estimated by leave-one-subject-out jackknife effects.

For the GLS consensus effect, bootstrap consensus draws are recomputed from the matched layer-effect draws and the fitted covariance structure. Percentile/basic consensus intervals are available.

When a bootstrap interval is selected in `run_omics_braid()`, analytic interval columns are preserved and the selected bootstrap interval is added explicitly. **Analytic SEs and p-values are not replaced.** This permits direct coverage comparison without changing the null-hypothesis machinery underneath the method.

## 9. Layer-level practical equivalence

For a user-defined smallest effect size of interest `Delta > 0`, each layer is tested for practical equivalence to `[-Delta, +Delta]` with TOST.

The local inferential states are:

- `equivalent`: practical equivalence is supported;
- `positive` / `negative`: a non-negligible directional effect is supported;
- `uncertain`: neither claim is established.

Multiplicity-adjusted states are also reported for screening across many entities, but local states are the default for describing an individual braid because adding unrelated entities should not alter that entity's geometric classification.

**A non-significant difference test is never treated as proof of equivalence.**

## 10. Cross-Omics GLS Trend Test

Version 0.2.0 retains the covariance-aware trajectory test introduced in v0.1.9 and does not alter it during confirmatory validation. The test replaces the earlier raw point-estimate slope rule with a covariance-aware trajectory test.

For an ordered omic sequence with layer index `L = 0, 1, ..., K-1`, effects are aligned to the dominant observed direction and modeled as

`y_aligned = beta_0 + beta_1 L + epsilon`,

with sampling covariance `V` and GLS estimator

`beta_hat = (X' V^-1 X)^-1 X' V^-1 y_aligned`.

The slope `beta_1` is interpreted as standardized-effect change per one-layer transition. A prespecified trajectory SESOI `Delta_T > 0` defines three inferential questions:

1. meaningful amplification: `beta_1 > +Delta_T`;
2. meaningful attenuation: `beta_1 < -Delta_T`;
3. practically flat trajectory: `-Delta_T < beta_1 < +Delta_T` via TOST.

Local and multiplicity-adjusted trajectory p-values/states are retained. Mixed-sign estimated braids are marked `not_applicable` for trajectory subtype interpretation and are handled by discordance/inversion logic. The default development value `Delta_T = 0.15` is not universal and must be sensitivity-tested and scientifically justified.

The trend test is descriptive across an ordered layer sequence. It does not prove temporal or causal propagation.

## 11. Hierarchical braid classification

Version 0.2.0 retains the v0.1.9 hierarchy separating **directional evidence**, **trajectory subtype**, **equivalence-dependent patterns**, and **lack of evidence**.

### Confirmed inversion

Supported positive and negative layer states coexist.

### Confirmed buffering

One or more supported directional upstream layers are followed only by downstream layers that **pass practical-equivalence testing**.

### Confirmed emergence

One or more upstream layers **pass practical-equivalence testing** before a supported downstream directional effect appears.

### Same-direction family

When all observed layers support the same direction:

- `attenuation` is assigned only when the GLS trend test confirms `beta_1 < -Delta_T`;
- `amplification` is assigned only when the GLS trend test confirms `beta_1 > +Delta_T`;
- `concordant_increase` / `concordant_decrease` is assigned when the trajectory is practically flat;
- if no trajectory subtype is confirmed, the broader concordant direction is retained with `pattern_status = direction_confirmed` rather than forcing a false attenuation/amplification call.

### Null-like outcomes

- `null_equivalent`: all observed layers pass practical-equivalence testing;
- `no_detectable_effect`: the joint-null omnibus test is not rejected, but practical equivalence has not been established;
- `uncertain`: available evidence does not justify a confirmatory braid label;
- `insufficient`: fewer than two usable layers.

This separation is fundamental: **absence of detected evidence and evidence of absence are not the same claim.**

## 12. Confirmatory versus suggestive geometry

Every entity receives a primary inferential `pattern` and `pattern_status`. OmicsBraid also reports `suggestive_pattern`, derived only from estimated effect magnitudes relative to the equivalence and trajectory margins.

A suggestive pattern is intended for cases where geometry is biologically interesting but precision is insufficient for confirmatory equivalence/trend claims. It must not be presented as a statistically confirmed pattern.

## 13. Geometric uncertainty propagation

For each entity, OmicsBraid draws from

`N(delta_hat, V_hat)`

and reclassifies each draw by practical effect regions and trajectory geometry. Frequencies quantify uncertainty in the **estimated braid geometry**. They are not Bayesian posterior probabilities and do not replace the confirmatory inferential label.

## 14. Pathway mode

RNA, protein, and metabolite measurements should not generally be treated as the same molecular entity. For heterogeneous molecular layers, pathway/reaction-level activity is usually the coherent unit. The built-in `score_pathways()` is a convenience baseline; externally validated pathway/activity scores may be supplied.

## 15. Validation design

The v0.2.0 confirmatory program freezes the v0.1.9 inferential/classification rules and concentrates simulation effort on calibration rather than additional classifier tuning.

The principal grid uses `n/group = 20, 40, 80, 160, 320`, normal versus heavy-tailed residuals, moderate cross-omic correlation (`rho = 0.30`), no imposed missingness, and at least 500 Monte-Carlo replicates per n-by-distribution cell. Each simulated dataset contains known joint-null, equal-effect, concordant, attenuation, amplification, inversion, buffering, and emergence entities.

The confirmatory targets are:

1. omnibus Type-I error under the joint null;
2. `Q_omics` Type-I error under equal nonzero effects;
3. matched-covariance versus false-independence calibration;
4. analytic, percentile, basic, and targeted BCa layer-CI coverage;
5. analytic versus percentile GLS-consensus CI coverage under equal effects;
6. attenuation/amplification trend power curves;
7. false meaningful-trend calls under equal effects;
8. equivalence power for truly zero buffering/emergence layers;
9. false equivalence for non-negligible layers;
10. decisive braid-classification safety;
11. scenario-execution reliability;
12. semi-synthetic recovery in real covariance/noise backgrounds after the simulation method is frozen;
13. external biological reproducibility after semi-synthetic validation.

Monte-Carlo confidence intervals are reported for proportions. Calibration is judged against compatibility with the nominal target and the full operating curve, not by demanding that a finite simulation point estimate equal exactly 0.05 or 0.95.

## 16. Scope restriction

Version 0.2.0 implements independent two-group contrasts. “Matched” refers to the same biological subjects being measured across omic layers; it does **not** mean paired case-control, crossover, longitudinal, or repeated-measures contrasts. Those designs require different effect estimators and resampling units.

## 17. Robust empirical calibration (v0.2.2)

The v0.2.1 confirmatory study isolated two remaining departures from nominal behavior under heavy-tailed simulated data: mild undercoverage of layer-level analytic intervals and modest anti-conservatism of the asymptotic chi-square reference for `Q_omics` in some sample-size cells. Version 0.2.2 therefore adds **optional empirical calibration without changing the validated braid classifier**.

### 17.1 Subject-level permutation omnibus

For the global null that every oriented omic-layer effect is zero, group labels may be permuted when the two-group design is exchangeable. A subject is the resampling unit; all available omic measurements from that subject remain linked.

For each layer `k`, the observed standardized effect is studentized by its analytic standard error, `z_k = g_k / SE(g_k)`. For every permuted label assignment the Hedges-g effect and its standard error are recomputed. Let `S_perm` be the covariance matrix of the jointly permuted studentized effect vectors and `W_perm = S_perm^{-1}` after positive-semidefinite regularization. The covariance-aware empirical omnibus statistic is

`T_omnibus = (z - mean_perm)^T W_perm (z - mean_perm)`.

The same quadratic form is evaluated for each centered joint permutation draw. The empirical p-value uses the finite-resampling correction

`p = (1 + #{T_b >= T_obs}) / (B + 1)`.

Because whole biological subjects are permuted, the joint permutation distribution retains cross-layer dependence directly. The test therefore avoids both an independence assumption and a chi-square reference approximation.

### 17.2 Centered-bootstrap omnibus

As a nonparametric comparator, a stratified matched-subject bootstrap is generated under the observed two-group design. Layer-specific bootstrap effects are centered by their bootstrap means. Their empirical covariance matrix defines a Mahalanobis/Wald metric; both the observed effect vector and centered bootstrap-null vectors are evaluated with that same covariance-aware quadratic form.

### 17.3 Centered-bootstrap heterogeneity test

The heterogeneity null is not the global zero null: it permits an unknown common non-zero standardized effect. Ordinary group-label permutation would therefore target the wrong null hypothesis.

Let `g_hat` be the observed K-layer effect vector and `g_b*` the matched-subject bootstrap vectors. For the centered-bootstrap comparator, each bootstrap layer is centered by its bootstrap mean. For either bootstrap null, the empirical covariance matrix `V*` defines `W* = (V*)^{-1}`. The fitted common effect for a vector `g` is the GLS estimate `theta = (1^T W* g)/(1^T W* 1)`, and the covariance-aware heterogeneity statistic is

`T_het = (g - theta 1)^T W* (g - theta 1)`.

The same statistic is evaluated for every empirical-null bootstrap draw, and the p-value uses `(1 + exceedances)/(B + 1)`. This is an empirical analogue of the GLS residual structure underlying `Q_omics`, while remaining separately reported so asymptotic and robust inference can be audited side-by-side.

### 17.4 Raw-data null-shift bootstrap for heterogeneity

Version 0.2.2 also implements a stronger heterogeneity calibration than the effect-level centered bootstrap. For each entity, OmicsBraid estimates the common oriented standardized effect across layers. Within each omic, the comparison-group observations are shifted by a constant so that the layer's Hedges-g effect equals that common target while leaving within-group residuals, heavy tails, and missing values unchanged. The adjusted multi-omic dataset therefore satisfies the fitted common-effect null by construction. Biological subjects are then resampled within groups as matched units. The covariance of those null bootstrap effects defines a GLS metric, and the covariance-aware residual quadratic heterogeneity statistic is recalculated for the observed vector and every null draw. This is the recommended empirical heterogeneity calibration in the v0.2.2 validation build.

### 17.5 Primary versus reported empirical inference

`run_omics_braid()` defaults to `empirical_tests = FALSE`. When empirical calibration is requested, asymptotic p-values remain the primary values unless the analyst explicitly sets `empirical_use_as_primary = TRUE`. Even then, the original asymptotic p-values are preserved in dedicated columns. This prevents silent changes in inferential method.

### 17.6 Exchangeability restriction

The subject-label permutation mode is valid only when the two-group labels are exchangeable under the null. Paired, longitudinal, clustered, covariate-adjusted, blocked, family-based, or repeated-measures designs require design-specific restricted permutations or alternative resampling schemes and are outside v0.2.2.

## 18. Final robust-calibration validation design

The v0.2.2 targeted validation freezes all classifier rules and tests only the new robust components. The primary grid uses `n/group = 20, 40, 80, 160, 320`, normal versus heavy-tailed t residuals, `rho = 0.30`, 500 outer Monte-Carlo replicates per cell, 300 matched covariance-bootstrap draws, 199 empirical calibration resamples, and targeted BCa coverage in the first 100 replicates per cell.

The simulated dataset contains three prespecified entities: a joint null `(0,0,0)`, an equal non-zero effect `(0.8,0.8,0.8)`, and a strong inversion `(1.0,0.8,-1.0)`. The benchmark compares asymptotic versus empirical Type-I error, inversion power, and analytic/percentile/basic/BCa interval coverage. High-frequency checkpoints are written to persistent internal-disk storage and only final results are synchronized to external project storage.
