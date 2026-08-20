# OmicsBraid input/output specification (v0.2.0)

## 1. Sample-level input

### Required

**Metadata**: one row per independent biological sample.

Minimum columns:

| column | meaning |
|---|---|
| `sample_id` | unique biological sample identifier |
| grouping variable | binary reference/comparison group |

**Assays**: named numeric matrices or data frames, one per omic.

- rows = features or precomputed pathway/activity scores;
- columns = sample IDs;
- row names = entity IDs;
- sample overlap across assays may be incomplete;
- assay values may contain `NA`;
- raw counts/raw MS intensities are not expected.

Examples: `RNA`, `Protein`, `Phosphoprotein`, `Metabolite`, `Lipid`.

### Optional

- feature annotation table;
- pathway mapping (`omic`, `feature_id`, `pathway`);
- precomputed pathway/activity matrices;
- prespecified orientation multipliers for omics whose biological direction has an inverse interpretation.

## 2. Summary-statistics input

Minimum columns:

| column | meaning |
|---|---|
| `entity` | common biological entity/pathway |
| `omic` | omic layer |
| `effect` | standardized effect on a common interpretation scale |
| `se` | standard error of the standardized effect |

Optional but recommended: `p_value`, `p_adj`, `conf_low`, `conf_high`, sample sizes, and entity-specific covariance matrices.

Do not directly combine incomparable coefficients such as RNA log2FC, methylation delta-beta, and arbitrary metabolite fold-change. Use a common standardized scale or sample-level mode.

## 3. Key user settings

- `omic_order`: scientific/display order of layers;
- `equivalence_margin`: SESOI for calling a layer practically negligible;
- `trajectory_margin`: SESOI for meaningful attenuation/amplification per one-layer transition;
- `alpha`: local inferential significance level;
- `bootstrap_B`: matched-subject bootstrap replicates;
- `bootstrap_shrinkage`: shrinkage of estimated cross-omic correlations;
- `ci_method`: `analytic`, `percentile`, or `bca` layer-effect interval;
- `integrated_ci_method`: `analytic`, `percentile`, or `basic` GLS-consensus interval;
- `ci_conf_level`: requested interval confidence level;
- `ci_min_boot`: minimum usable bootstrap draws before a bootstrap interval replaces the analytic display interval;
- `orientation`: prespecified +1/-1 direction harmonization.

Neither SESOI has a universal value; both require scientific justification and sensitivity analysis.

## 4. Primary output tables

### `layer_effects.csv`

One row per entity x omic with Hedges' g, analytic SE/p value, group sample sizes, and the selected display CI. When bootstrap intervals are requested, analytic CI columns are preserved alongside bootstrap CI columns and the number of usable bootstrap draws.

### `integrated_effects.csv`

One row per entity containing:

- GLS consensus standardized effect and CI;
- local/adjusted consensus p values;
- multivariate `W_omnibus` and omnibus p values;
- `Q_omics` and heterogeneity p values;
- descriptive `I2_omics`;
- weighted directional agreement;
- evidence-qualified directional agreement;
- GLS-weight and covariance-condition diagnostics;
- covariance mode;
- analytic versus bootstrap consensus interval columns when a robust consensus CI is requested.


### `effect_intervals` / `consensus_intervals`

When bootstrap intervals are requested, the result object also stores interval-diagnostic tables containing bootstrap mean, SD, bias, number of usable draws, selected limits, and BCa bias-correction/acceleration diagnostics where applicable.

### `equivalence_results.csv`

Adds:

- omic-specific equivalence margin;
- lower/upper TOST p values;
- local and adjusted TOST evidence;
- local and adjusted practical states: `positive`, `negative`, `equivalent`, `uncertain`.

### `trend_results.csv`

One row per entity containing the ordered covariance-aware GLS trajectory:

- aligned direction;
- slope and SE;
- 95% CI;
- p value for zero slope;
- one-sided p values for meaningful attenuation/amplification relative to `trajectory_margin`;
- TOST p value for a practically flat trajectory;
- local and multiplicity-adjusted trend state (`attenuation`, `amplification`, `flat`, `uncertain`, or `not_applicable` for mixed-sign braids);
- covariance mode used for the trend.

### `braid_classification.csv`

One row per entity containing:

- primary inferential `pattern`;
- `pattern_status`: `confirmed`, `direction_confirmed`, `no_evidence`, `unresolved`, or `insufficient`;
- `suggestive_pattern` based on effect geometry only;
- human-readable `interpretation_label`;
- confirmatory basis;
- counts of positive/negative/equivalent/uncertain layer states;
- local omnibus evidence;
- trend diagnostics;
- geometric Monte Carlo support/entropy diagnostics.

Primary pattern labels include:

- `concordant_increase`, `concordant_decrease`;
- `attenuation`, `amplification`;
- `buffering`, `emergence`;
- `inversion`;
- `null_equivalent`;
- `no_detectable_effect`;
- `uncertain`, `insufficient`.

### `pattern_probabilities.csv`

Long table of uncertainty-propagation frequencies for practical **geometric** patterns. These frequencies are not Bayesian posterior probabilities.

### `covariance_diagnostics.csv`

For covariance-aware sample-level analyses:

- number of omics;
- minimum usable bootstrap pairs;
- maximum absolute estimated cross-omic correlation;
- magnitude of PSD repair.

## 5. Primary figures

- **Omics Evidence Forest**: layer effects + CIs + consensus effect + omnibus/Q/I2 summary.
- **Effect Braid**: ordered effects, CIs, inferential layer states, equivalence band, hierarchical pattern/status, and trend state.
- **Concordance Map**: omnibus evidence versus evidence-qualified directional agreement by default.
- **Braid Heatmap**: compact layer-specific standardized effects across entities.

## 6. Interpretation hierarchy

For an entity/pathway:

1. inspect layer effects and CIs;
2. inspect the multivariate omnibus test;
3. inspect `Q_omics` and heterogeneity;
4. inspect layer equivalence/difference states;
5. inspect the ordered GLS trend test;
6. inspect the primary braid pattern and `pattern_status`;
7. inspect suggestive geometry and uncertainty propagation;
8. interpret the GLS consensus last when heterogeneity is substantial.

A `no_detectable_effect` result is not equivalent to `null_equivalent`. A `suggestive_pattern` is not a confirmatory pattern.
