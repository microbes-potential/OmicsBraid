#' Plot an Omics Evidence Forest
#' @param result An `omics_braid_result`.
#' @param entity Entity/pathway to display.
#' @param omic_order Optional omic order.
#' @return A ggplot object.
#' @export
plot_evidence_forest <- function(result, entity, omic_order = result$settings$omic_order) {
  d <- result$effects[result$effects$entity == entity & result$effects$omic %in% omic_order, , drop = FALSE]
  if (!nrow(d)) .stopf("Entity '%s' not found in layer-specific effects.", entity)
  d <- d[match(intersect(omic_order, d$omic), d$omic), , drop = FALSE]
  d$label <- d$omic
  d$type <- "Layer"
  it <- result$integrated[result$integrated$entity == entity, , drop = FALSE]
  if (nrow(it)) {
    di <- data.frame(entity = entity, omic = "Integrated", effect = it$integrated_effect,
                     se = it$integrated_se, conf_low = it$conf_low, conf_high = it$conf_high,
                     label = "Integrated", type = "Integrated", stringsAsFactors = FALSE)
    d <- rbind(d[, intersect(names(d), names(di)), drop = FALSE], di[, intersect(names(d), names(di)), drop = FALSE])
  }
  d$label <- factor(d$label, levels = rev(c(intersect(omic_order, as.character(d$label)), if ("Integrated" %in% d$label) "Integrated" else character())))
  subtitle <- if (nrow(it)) sprintf("Omnibus p = %.3g; Qomics = %.2f (heterogeneity p = %.3g); I2-omics = %.1f%%",
                                   it$p_omnibus, it$Q_omics, it$p_heterogeneity, it$I2_omics) else NULL
  ggplot2::ggplot(d, ggplot2::aes(x = effect, y = label, xmin = conf_low, xmax = conf_high, shape = type)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_pointrange(size = 0.45) +
    ggplot2::scale_shape_manual(values = c(Layer = 16, Integrated = 18)) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = entity, subtitle = subtitle, x = "Standardized effect", y = NULL, shape = NULL)
}

#' Plot an Effect Braid
#' @param result An `omics_braid_result`.
#' @param entity Entity/pathway to display.
#' @param omic_order Optional omic order.
#' @param show_ci Show confidence intervals.
#' @param show_equivalence_region Shade the negligible-effect region.
#' @return A ggplot object.
#' @export
plot_effect_braid <- function(result, entity, omic_order = result$settings$omic_order,
                              show_ci = TRUE, show_equivalence_region = TRUE) {
  d <- result$equivalence[result$equivalence$entity == entity & result$equivalence$omic %in% omic_order, , drop = FALSE]
  if (!nrow(d)) .stopf("Entity '%s' not found.", entity)
  d <- d[match(intersect(omic_order, d$omic), d$omic), , drop = FALSE]
  d$omic <- factor(d$omic, levels = omic_order)
  cl <- result$classification[result$classification$entity == entity, , drop = FALSE]
  subtitle <- if (nrow(cl)) {
    label <- if ("interpretation_label" %in% names(cl) && nzchar(cl$interpretation_label[1])) cl$interpretation_label[1] else cl$pattern[1]
    status <- if ("pattern_status" %in% names(cl) && nzchar(cl$pattern_status[1])) paste0("; status ", cl$pattern_status[1]) else ""
    sup <- if ("deterministic_pattern_support" %in% names(cl) && is.finite(cl$deterministic_pattern_support[1])) {
      sprintf("; geometry support %.2f", cl$deterministic_pattern_support[1])
    } else ""
    tr <- if ("trend_state_local" %in% names(cl) && nzchar(cl$trend_state_local[1])) {
      paste0("; trend ", cl$trend_state_local[1])
    } else ""
    paste0("Pattern: ", label, status, tr, sup)
  } else NULL
  p <- ggplot2::ggplot(d, ggplot2::aes(x = omic, y = effect, group = 1))
  if (show_equivalence_region) {
    band <- data.frame(xmin = seq_len(nrow(d)) - 0.5, xmax = seq_len(nrow(d)) + 0.5,
                       ymin = -d$equiv_margin, ymax = d$equiv_margin)
    p <- p + ggplot2::geom_rect(data = band,
                                ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                                inherit.aes = FALSE, alpha = 0.07)
  }
  p <- p +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(shape = state), size = 3) +
    ggplot2::scale_shape_manual(values = c(positive = 24, negative = 25, equivalent = 21, uncertain = 4)) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = entity, subtitle = subtitle, x = NULL, y = "Standardized effect", shape = "State")
  if (show_ci && all(c("conf_low","conf_high") %in% names(d))) p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = conf_low, ymax = conf_high), width = 0.12)
  p
}

#' Plot concordance versus integrated significance
#' @param result An `omics_braid_result`.
#' @param label_top Number of most significant entities to label.
#' @param metric Concordance metric: evidence-qualified directional agreement (default), raw weighted directional agreement, or descriptive I2-based consistency.
#' @param significance Evidence axis: multivariate omnibus (default) or GLS consensus-effect significance.
#' @return A ggplot object.
#' @export
plot_concordance_map <- function(result, label_top = 0L,
                                 metric = c("evidence_direction_agreement", "direction_agreement", "i2_consistency"),
                                 significance = c("omnibus", "consensus")) {
  metric <- match.arg(metric)
  significance <- match.arg(significance)
  d <- result$integrated
  pcol <- if (significance == "omnibus") "p_omnibus_adj" else "p_adj"
  d$minus_log10_p <- -log10(pmax(d[[pcol]], .Machine$double.xmin))
  if (metric == "evidence_direction_agreement") {
    d$concordance <- d$evidence_qualified_direction_agreement
    ylab <- "Evidence-qualified directional agreement"
  } else if (metric == "direction_agreement") {
    d$concordance <- d$direction_agreement
    ylab <- "Raw weighted directional alignment"
  } else {
    d$concordance <- 1 - pmin(d$I2_omics, 100) / 100
    ylab <- "I2-based consistency (descriptive)"
  }
  p <- ggplot2::ggplot(d, ggplot2::aes(x = minus_log10_p, y = concordance, size = n_omics)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(x = if (significance == "omnibus") expression(-log[10](adjusted~p[omnibus])) else expression(-log[10](adjusted~p[consensus])),
                  y = ylab, size = "Omics")
  if (label_top > 0L) {
    ii <- utils::head(order(d[[pcol]]), label_top)
    p <- p + ggplot2::geom_text(data = d[ii, , drop = FALSE], ggplot2::aes(label = entity), nudge_y = 0.03, check_overlap = TRUE, size = 3)
  }
  p
}

#' Plot a braid heatmap across entities and omics
#' @param result An `omics_braid_result`.
#' @param entities Optional entities; default top 25 by integrated adjusted p-value.
#' @param omic_order Optional omic order.
#' @return A ggplot object.
#' @export
plot_braid_heatmap <- function(result, entities = NULL, omic_order = result$settings$omic_order) {
  if (is.null(entities)) entities <- utils::head(result$integrated$entity[order(result$integrated$p_omnibus_adj)], 25L)
  d <- result$effects[result$effects$entity %in% entities & result$effects$omic %in% omic_order, , drop = FALSE]
  d$entity <- factor(d$entity, levels = rev(entities))
  d$omic <- factor(d$omic, levels = omic_order)
  ggplot2::ggplot(d, ggplot2::aes(x = omic, y = entity, fill = effect)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(midpoint = 0, name = "Effect") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(x = NULL, y = NULL)
}
