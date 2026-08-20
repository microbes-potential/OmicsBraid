# Package index

## Complete workflows

- [`run_omics_braid()`](https://microbes-potential.github.io/OmicsBraid/reference/run_omics_braid.md)
  : Run the complete OmicsBraid workflow on sample-level data
- [`run_omics_braid_summary()`](https://microbes-potential.github.io/OmicsBraid/reference/run_omics_braid_summary.md)
  : Run OmicsBraid from externally estimated summary statistics
- [`braid_results_table()`](https://microbes-potential.github.io/OmicsBraid/reference/braid_results_table.md)
  : Create a one-row-per-entity master result table
- [`write_omics_braid()`](https://microbes-potential.github.io/OmicsBraid/reference/write_omics_braid.md)
  : Export OmicsBraid results

## Data input and harmonization

- [`omics_braid_data()`](https://microbes-potential.github.io/OmicsBraid/reference/omics_braid_data.md)
  : Construct an OmicsBraid data object
- [`as_omics_braid_data()`](https://microbes-potential.github.io/OmicsBraid/reference/as_omics_braid_data.md)
  : Coerce supported objects to OmicsBraid data
- [`validate_omics_braid_data()`](https://microbes-potential.github.io/OmicsBraid/reference/validate_omics_braid_data.md)
  : Validate an OmicsBraid data object
- [`read_omics_braid()`](https://microbes-potential.github.io/OmicsBraid/reference/read_omics_braid.md)
  : Read OmicsBraid input files
- [`harmonize_entities()`](https://microbes-potential.github.io/OmicsBraid/reference/harmonize_entities.md)
  : Harmonize assay-specific feature identifiers to common entities
- [`score_pathways()`](https://microbes-potential.github.io/OmicsBraid/reference/score_pathways.md)
  : Score pathways within each omic layer
- [`orient_omics()`](https://microbes-potential.github.io/OmicsBraid/reference/orient_omics.md)
  : Harmonize the sign orientation of omic-layer effects

## Effect estimation and covariance

- [`estimate_effects()`](https://microbes-potential.github.io/OmicsBraid/reference/estimate_effects.md)
  : Estimate layer-specific standardized effects
- [`bootstrap_effect_covariance()`](https://microbes-potential.github.io/OmicsBraid/reference/bootstrap_effect_covariance.md)
  : Estimate cross-omic effect covariance by matched-subject bootstrap
- [`bootstrap_effect_intervals()`](https://microbes-potential.github.io/OmicsBraid/reference/bootstrap_effect_intervals.md)
  : Bootstrap confidence intervals for layer-specific standardized
  effects
- [`bootstrap_consensus_intervals()`](https://microbes-potential.github.io/OmicsBraid/reference/bootstrap_consensus_intervals.md)
  : Bootstrap confidence intervals for GLS consensus effects

## Cross-omic inference

- [`integrate_effects()`](https://microbes-potential.github.io/OmicsBraid/reference/integrate_effects.md)
  : Integrate effects across omic layers and quantify heterogeneity
- [`empirical_omics_tests()`](https://microbes-potential.github.io/OmicsBraid/reference/empirical_omics_tests.md)
  : Empirically calibrate OmicsBraid omnibus and heterogeneity tests
- [`test_equivalence()`](https://microbes-potential.github.io/OmicsBraid/reference/test_equivalence.md)
  : Test practical equivalence to a negligible effect region
- [`test_braid_trend()`](https://microbes-potential.github.io/OmicsBraid/reference/test_braid_trend.md)
  : Test ordered cross-omic effect trajectories with generalized least
  squares
- [`classify_braids()`](https://microbes-potential.github.io/OmicsBraid/reference/classify_braids.md)
  : Classify cross-omic braid patterns using inferential evidence
- [`braid_pattern_probabilities()`](https://microbes-potential.github.io/OmicsBraid/reference/braid_pattern_probabilities.md)
  : Quantify uncertainty in effect-braid geometry

## Visualization

- [`plot_evidence_forest()`](https://microbes-potential.github.io/OmicsBraid/reference/plot_evidence_forest.md)
  : Plot an Omics Evidence Forest
- [`plot_effect_braid()`](https://microbes-potential.github.io/OmicsBraid/reference/plot_effect_braid.md)
  : Plot an Effect Braid
- [`plot_concordance_map()`](https://microbes-potential.github.io/OmicsBraid/reference/plot_concordance_map.md)
  : Plot concordance versus integrated significance
- [`plot_braid_heatmap()`](https://microbes-potential.github.io/OmicsBraid/reference/plot_braid_heatmap.md)
  : Plot a braid heatmap across entities and omics

## Simulation

- [`simulate_braid_data()`](https://microbes-potential.github.io/OmicsBraid/reference/simulate_braid_data.md)
  : Simulate multi-omics data with known braid patterns
