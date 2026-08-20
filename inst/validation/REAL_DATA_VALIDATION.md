# Real-data validation plan for OmicsBraid

Real cohorts are **biological validation**, not statistical ground truth. They should be used only after null/alternative simulations and semi-synthetic injection tests establish operating characteristics.

## Primary development cohort: CPTAC AML multi-omic atlas (2026)

Recommended because the published cohort contains 173 treatment-naive participants spanning genomics, methylomics, transcriptomics, proteomics, PTMs, metabolomics and lipidomics. The article reports RNA-seq for 172 tumors, proteomics including PTMs for 162 tumors, metabolomics for 91 tumors and lipidomics for 96 tumors (with 97 participants described as having metabolomics and/or lipidomics in the cohort overview), creating both a rich matched subset and realistic incomplete-modality structure.

Suggested OmicsBraid use:

1. Define a biologically defensible binary contrast (for example, a pre-specified molecular subtype or mutation-defined group with adequate sample size).
2. Use analysis-ready RNA, protein, PTM and metabolite/lipid pathway scores.
3. Run entity-matched RNA -> protein -> PTM analyses where IDs are truly linked.
4. Run pathway-level RNA -> protein -> metabolite analyses for cross-chemistry integration.
5. Compare independence-assumption results with matched-subject bootstrap covariance results.
6. Test whether pathways with high `Q_omics` correspond to interpretable cross-layer discordance rather than technical artifacts.
7. Perform split-sample stability analysis without changing the locked braid definitions.

Data resources include NCI Proteomic Data Commons studies beginning with PDC000554 (proteome) and associated AML studies through the lipidome study PDC000562; genomic/transcriptomic components are linked through CPTAC/GDC resources.

Primary article: Chu SCA et al. *Integrated proteogenomic and metabolomic profiling of acute myeloid leukemia*. Nature Cancer (2026), article s43018-026-01175-6.

## Secondary validation: CPTAC glioblastoma

The CPTAC GBM discovery resource provides a distinct solid-tumor setting with RNA/proteome/phosphoproteome measurements and published multi-omic analyses. It is useful for demonstrating that OmicsBraid behavior is not specific to a hematologic malignancy.

Suggested use:

- entity-level RNA -> protein -> phosphoprotein evidence forests;
- pathway-level analyses where metabolite measurements are available in the chosen release/companion resource;
- sensitivity of `Q_omics` and braid classifications to sample composition.

## Replication-oriented validation

For a final methods manuscript, add a disease with both discovery and independent confirmatory cohorts. Lock the effect definitions, equivalence margin strategy and braid rules in the discovery cohort, then test whether high-confidence braid patterns reproduce in the confirmatory cohort.

## Minimum real-data reporting

For every cohort report:

- exact sample overlap across omics;
- preprocessing performed before OmicsBraid;
- number of entities/pathways tested;
- equivalence margin and its scientific justification;
- bootstrap B, correlation shrinkage and minimum pair count;
- fraction of entities whose covariance required PSD repair;
- integrated p-value/FDR distribution;
- `Q_omics` calibration checks/sensitivity analysis;
- pattern counts and Monte Carlo pattern probabilities;
- results under independence assumption as a sensitivity comparator;
- stability under subject bootstrap/split samples;
- whether conclusions persist when alternative valid pathway-scoring methods are used.
