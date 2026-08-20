#' Create a one-row-per-entity master result table
#'
#' @param result An `omics_braid_result`.
#' @return Data frame merging integrated inference, braid classification, and compact layer summaries.
#' @export
braid_results_table <- function(result) {
  x <- merge(result$integrated, result$classification, by = "entity", all = TRUE, sort = FALSE)
  state_summary <- lapply(split(result$equivalence, result$equivalence$entity), function(d) {
    data.frame(entity = d$entity[1],
               layer_states = paste(paste0(d$omic, ":", d$state), collapse = ";"),
               layer_effects = paste(paste0(d$omic, ":", sprintf("%.3f", d$effect)), collapse = ";"),
               stringsAsFactors = FALSE)
  })
  ss <- if (length(state_summary)) do.call(rbind, state_summary) else data.frame(entity = character())
  merge(x, ss, by = "entity", all.x = TRUE, sort = FALSE)
}

#' Export OmicsBraid results
#'
#' @param result An `omics_braid_result`.
#' @param dir Output directory.
#' @param save_plots Save overview PDF figures.
#' @param top_n Number of top entities for heatmap and individual plots.
#' @return Invisibly returns the output directory.
#' @export
write_omics_braid <- function(result, dir, save_plots = TRUE, top_n = 20L) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(result$effects, file.path(dir, "layer_effects.csv"), row.names = FALSE)
  utils::write.csv(result$integrated, file.path(dir, "integrated_effects.csv"), row.names = FALSE)
  if (!is.null(result$empirical_tests)) utils::write.csv(result$empirical_tests, file.path(dir, "empirical_test_calibration.csv"), row.names = FALSE)
  if (!is.null(result$effect_intervals)) utils::write.csv(result$effect_intervals, file.path(dir, "layer_bootstrap_intervals.csv"), row.names = FALSE)
  if (!is.null(result$consensus_intervals)) utils::write.csv(result$consensus_intervals, file.path(dir, "consensus_bootstrap_intervals.csv"), row.names = FALSE)
  utils::write.csv(result$equivalence, file.path(dir, "equivalence_results.csv"), row.names = FALSE)
  if (!is.null(result$trend)) utils::write.csv(result$trend, file.path(dir, "trend_results.csv"), row.names = FALSE)
  utils::write.csv(result$classification, file.path(dir, "braid_classification.csv"), row.names = FALSE)
  utils::write.csv(result$pattern_probabilities, file.path(dir, "pattern_probabilities.csv"), row.names = FALSE)
  utils::write.csv(braid_results_table(result), file.path(dir, "master_results.csv"), row.names = FALSE)
  saveRDS(result$settings, file.path(dir, "analysis_settings.rds"))
  utils::capture.output(utils::str(result$settings), file = file.path(dir, "analysis_manifest.txt"))
  if (!is.null(result$covariance)) {
    saveRDS(result$covariance, file.path(dir, "cross_omic_covariance.rds"))
    if (!is.null(result$covariance$diagnostics)) utils::write.csv(result$covariance$diagnostics, file.path(dir, "covariance_diagnostics.csv"), row.names = FALSE)
  }
  ps <- if (!is.null(result$data)) attr(result$data, "pathway_scoring") else NULL
  if (!is.null(ps) && !is.null(ps$mapping_coverage)) utils::write.csv(ps$mapping_coverage, file.path(dir, "pathway_mapping_coverage.csv"), row.names = FALSE)
  if (save_plots) {
    pdir <- file.path(dir, "plots"); dir.create(pdir, showWarnings = FALSE)
    ggplot2::ggsave(file.path(pdir, "concordance_map.pdf"), plot_concordance_map(result, label_top = min(10L, top_n)), width = 7, height = 5)
    top <- utils::head(result$integrated$entity[order(result$integrated$p_omnibus_adj)], top_n)
    ggplot2::ggsave(file.path(pdir, "braid_heatmap.pdf"), plot_braid_heatmap(result, top), width = 7, height = max(5, 0.25 * length(top) + 2))
    for (e in utils::head(top, min(10L, length(top)))) {
      safe <- gsub("[^A-Za-z0-9_.-]", "_", e)
      ggplot2::ggsave(file.path(pdir, paste0(safe, "_forest.pdf")), plot_evidence_forest(result, e), width = 7, height = 4.5)
      ggplot2::ggsave(file.path(pdir, paste0(safe, "_braid.pdf")), plot_effect_braid(result, e), width = 7, height = 4.5)
    }
  }
  invisible(normalizePath(dir, mustWork = FALSE))
}
