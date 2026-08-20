#' Simulate multi-omics data with known braid patterns
#'
#' Generates matched sample-level data for method development and validation.
#' Residuals are correlated across omic layers for each entity, allowing the
#' covariance, heterogeneity, equivalence, and classification procedures to be
#' tested against known cross-layer effects.
#'
#' @param n_per_group Default samples per group when `n_reference` and `n_comparison` are not supplied.
#' @param n_reference Optional reference-group sample size.
#' @param n_comparison Optional comparison-group sample size.
#' @param omics Ordered omic names.
#' @param patterns Named list of true standardized mean shifts, one numeric vector per entity.
#' @param rho Scalar equicorrelation or an omic-by-omic residual correlation matrix.
#' @param missing_rate Independent value-level missingness probability.
#' @param modality_missing_rate Probability of an entire subject modality being absent; scalar or named by omic.
#' @param residual_distribution `"normal"` or heavy-tailed `"t"` residuals.
#' @param t_df Degrees of freedom for t residuals; must exceed 2 for finite variance.
#' @param seed Random seed.
#' @return List with `data` (`omics_braid_data`) and `truth` table.
#' @export
simulate_braid_data <- function(n_per_group = 60L, n_reference = NULL, n_comparison = NULL,
                                omics = c("RNA", "Protein", "Metabolite"),
                                patterns = NULL, rho = 0.4,
                                missing_rate = 0, modality_missing_rate = 0,
                                residual_distribution = c("normal", "t"), t_df = 5,
                                seed = 1L) {
  residual_distribution <- match.arg(residual_distribution)
  if (is.null(n_reference)) n_reference <- n_per_group
  if (is.null(n_comparison)) n_comparison <- n_per_group
  n_reference <- as.integer(n_reference); n_comparison <- as.integer(n_comparison)
  if (n_reference < 3L || n_comparison < 3L) .stopf("Both simulated groups require at least 3 samples.")
  if (!length(omics) || anyDuplicated(omics)) .stopf("'omics' must contain unique layer names.")
  K <- length(omics)
  if (K < 2L) .stopf("Simulation requires at least two omic layers.")

  if (is.null(patterns)) {
    patterns <- list(
      concordant_up = seq(1.0, 0.8, length.out = K),
      attenuation = seq(1.4, 0.5, length.out = K),
      amplification = seq(0.5, 1.4, length.out = K),
      inversion = c(rep(0.9, K - 1L), -1.0),
      buffering = c(1.1, rep(0, K - 1L)),
      emergence = c(rep(0, K - 1L), 1.1),
      null = rep(0, K)
    )
  }
  if (any(lengths(patterns) != K)) .stopf("Every pattern vector must have one value per omic.")
  patterns <- lapply(patterns, function(v) {
    if (!is.null(names(v))) v <- v[omics]
    v <- as.numeric(v)
    if (any(!is.finite(v))) .stopf("Pattern vectors must contain finite effects for every omic.")
    v
  })

  if (length(rho) == 1L) {
    if (!is.finite(rho) || rho <= -1/(K-1) || rho >= 1) .stopf("Scalar 'rho' must define a positive-definite equicorrelation matrix.")
    Sigma <- matrix(rho, K, K); diag(Sigma) <- 1
  } else {
    Sigma <- as.matrix(rho)
    if (!all(dim(Sigma) == c(K,K))) .stopf("Matrix 'rho' must be K x K where K = length(omics).")
    if (any(abs(diag(Sigma) - 1) > 1e-8)) .stopf("A supplied residual correlation matrix must have unit diagonal.")
    if (max(abs(Sigma - t(Sigma))) > 1e-8) .stopf("A supplied residual correlation matrix must be symmetric.")
    if (min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values) <= 0) .stopf("A supplied residual correlation matrix must be positive definite.")
  }

  if (missing_rate < 0 || missing_rate >= 1) .stopf("'missing_rate' must be in [0,1).")
  mm <- if (length(modality_missing_rate) == 1L) rep(as.numeric(modality_missing_rate), K) else {
    if (is.null(names(modality_missing_rate)) || !all(omics %in% names(modality_missing_rate))) .stopf("Named 'modality_missing_rate' must cover every omic.")
    as.numeric(modality_missing_rate[omics])
  }
  if (any(mm < 0 | mm >= 1)) .stopf("All modality missingness probabilities must be in [0,1).")
  if (residual_distribution == "t" && (!is.finite(t_df) || t_df <= 2)) .stopf("'t_df' must exceed 2.")

  E <- length(patterns)
  n <- n_reference + n_comparison
  ids <- sprintf("S%03d", seq_len(n))
  group <- c(rep("Control", n_reference), rep("Disease", n_comparison))
  assays <- stats::setNames(lapply(omics, function(o) matrix(NA_real_, E, n, dimnames = list(names(patterns), ids))), omics)
  set.seed(seed)
  for (e in seq_along(patterns)) {
    Z <- .rmvn_psd(n, rep(0, K), Sigma)
    if (residual_distribution == "t") {
      scale_mix <- sqrt(stats::rchisq(n, df = t_df) / t_df)
      Z <- Z / scale_mix
      Z <- Z / sqrt(t_df / (t_df - 2))
    }
    shift <- ifelse(group == "Disease", 1, 0)
    for (k in seq_len(K)) assays[[k]][e, ] <- Z[,k] + shift * patterns[[e]][k]
  }

  if (missing_rate > 0) {
    for (k in seq_along(assays)) {
      miss <- matrix(stats::runif(length(assays[[k]])) < missing_rate, nrow = nrow(assays[[k]]))
      assays[[k]][miss] <- NA_real_
    }
  }
  if (any(mm > 0)) {
    for (k in seq_along(assays)) {
      drop <- stats::runif(ncol(assays[[k]])) < mm[k]
      # Retain enough subjects to keep the simulated object meaningful.
      if (sum(!drop & group == "Control") >= 3L && sum(!drop & group == "Disease") >= 3L) {
        assays[[k]] <- assays[[k]][, !drop, drop = FALSE]
      }
    }
  }

  meta <- data.frame(sample_id = ids, group = group, stringsAsFactors = FALSE)
  truth <- do.call(rbind, lapply(names(patterns), function(e) {
    data.frame(entity = e, omic = omics, true_effect = patterns[[e]], stringsAsFactors = FALSE)
  }))
  list(data = omics_braid_data(assays, meta), truth = truth,
       simulation_settings = list(n_reference = n_reference, n_comparison = n_comparison,
                                  omics = omics, rho = Sigma, missing_rate = missing_rate,
                                  modality_missing_rate = stats::setNames(mm, omics),
                                  residual_distribution = residual_distribution, t_df = t_df, seed = seed))
}
