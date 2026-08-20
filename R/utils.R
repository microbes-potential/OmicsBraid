`%||%` <- function(x, y) if (is.null(x)) y else x

.stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)
.warnf <- function(fmt, ...) warning(sprintf(fmt, ...), call. = FALSE)

.is_scalar_string <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)

.as_numeric_matrix <- function(x, name = "assay") {
  if (is.data.frame(x)) {
    default_rn <- identical(rownames(x), as.character(seq_len(nrow(x))))
    if (ncol(x) > 1L && !is.numeric(x[[1L]]) && default_rn) {
      rn <- as.character(x[[1L]])
      x <- x[-1L]
      rownames(x) <- rn
    }
    x <- as.matrix(x)
  }
  if (!is.matrix(x)) .stopf("%s must be a matrix or data.frame.", name)
  storage.mode(x) <- "double"
  if (any(is.infinite(x), na.rm = TRUE)) .stopf("%s contains infinite values. Replace or resolve +/-Inf during assay-specific preprocessing before OmicsBraid.", name)
  if (is.null(rownames(x)) || any(!nzchar(rownames(x)))) .stopf("%s must have non-empty feature row names.", name)
  if (is.null(colnames(x)) || any(!nzchar(colnames(x)))) .stopf("%s must have non-empty sample column names.", name)
  if (anyDuplicated(rownames(x))) .stopf("%s contains duplicated feature IDs.", name)
  if (anyDuplicated(colnames(x))) .stopf("%s contains duplicated sample IDs.", name)
  x
}

.hedges_J <- function(df) {
  out <- rep(NA_real_, length(df))
  ok <- is.finite(df) & df > 1
  out[ok] <- exp(lgamma(df[ok] / 2) - 0.5 * log(df[ok] / 2) - lgamma((df[ok] - 1) / 2))
  out
}

.row_mean_na <- function(x) {
  n <- rowSums(!is.na(x))
  s <- rowSums(x, na.rm = TRUE)
  ans <- s / n
  ans[n == 0L] <- NA_real_
  ans
}

.row_var_na <- function(x) {
  n <- rowSums(!is.na(x))
  m <- .row_mean_na(x)
  z <- x - m
  z[is.na(z)] <- 0
  ss <- rowSums(z * z)
  ans <- ss / (n - 1)
  ans[n < 2L] <- NA_real_
  ans
}

.hedges_g_matrix <- function(mat, group, reference, comparison, min_n = 3L) {
  i0 <- which(group == reference)
  i1 <- which(group == comparison)
  if (!length(i0) || !length(i1)) return(rep(NA_real_, nrow(mat)))
  x0 <- mat[, i0, drop = FALSE]
  x1 <- mat[, i1, drop = FALSE]
  n0 <- rowSums(!is.na(x0)); n1 <- rowSums(!is.na(x1))
  m0 <- .row_mean_na(x0); m1 <- .row_mean_na(x1)
  v0 <- .row_var_na(x0); v1 <- .row_var_na(x1)
  df <- n0 + n1 - 2
  sp2 <- ((n0 - 1) * v0 + (n1 - 1) * v1) / df
  d <- (m1 - m0) / sqrt(sp2)
  g <- .hedges_J(df) * d
  g[n0 < min_n | n1 < min_n | !is.finite(g) | sp2 <= 0] <- NA_real_
  g
}

.make_psd <- function(M, eps = 1e-8) {
  M <- (M + t(M)) / 2
  if (nrow(M) == 1L) return(matrix(max(M[1,1], eps), 1, 1, dimnames = dimnames(M)))
  ee <- eigen(M, symmetric = TRUE)
  floorv <- max(eps, max(abs(ee$values), na.rm = TRUE) * eps)
  vals <- pmax(ee$values, floorv)
  out <- ee$vectors %*% (vals * t(ee$vectors))
  out <- (out + t(out)) / 2
  dimnames(out) <- dimnames(M)
  out
}

.safe_inverse <- function(M, tol = 1e-10) {
  M <- .make_psd(M, eps = tol)
  ee <- eigen(M, symmetric = TRUE)
  cutoff <- max(ee$values) * tol
  invv <- ifelse(ee$values > cutoff, 1 / ee$values, 0)
  ee$vectors %*% (invv * t(ee$vectors))
}

.rmvn_psd <- function(n, mu, Sigma) {
  Sigma <- .make_psd(Sigma)
  ee <- eigen(Sigma, symmetric = TRUE)
  A <- ee$vectors %*% diag(sqrt(pmax(ee$values, 0)), nrow = length(ee$values))
  Z <- matrix(stats::rnorm(n * length(mu)), nrow = n)
  sweep(Z %*% t(A), 2L, mu, "+")
}


.validate_effect_table <- function(effects, require_p = FALSE) {
  if (!is.data.frame(effects)) .stopf("'effects' must be a data.frame.")
  req <- c("entity", "omic", "effect", "se")
  if (isTRUE(require_p)) req <- c(req, "p_value")
  if (!all(req %in% names(effects))) .stopf("effects must contain columns: %s", paste(req, collapse = ", "))
  effects$entity <- as.character(effects$entity)
  effects$omic <- as.character(effects$omic)
  if (anyNA(effects$entity) || any(!nzchar(effects$entity)) ||
      anyNA(effects$omic) || any(!nzchar(effects$omic))) {
    .stopf("'entity' and 'omic' must be non-missing, non-empty strings.")
  }
  key <- paste(effects$entity, effects$omic, sep = "\r")
  if (anyDuplicated(key)) {
    bad <- unique(key[duplicated(key) | duplicated(key, fromLast = TRUE)])
    shown <- gsub("\r", " / ", utils::head(bad, 5L), fixed = TRUE)
    .stopf("Each entity x omic combination must occur once. Duplicates include: %s", paste(shown, collapse = ", "))
  }
  if (!is.numeric(effects$effect) || !is.numeric(effects$se)) .stopf("'effect' and 'se' must be numeric.")
  bad_se <- is.finite(effects$effect) & (!is.finite(effects$se) | effects$se <= 0)
  if (any(bad_se)) .stopf("Every finite effect must have a finite standard error > 0.")
  invisible(TRUE)
}

.resolve_margin <- function(margin, omics) {
  if (length(margin) == 1L) return(rep(as.numeric(margin), length(omics)))
  if (is.null(names(margin))) .stopf("When 'margin' has length > 1, it must be named by omic.")
  if (!all(omics %in% names(margin))) .stopf("No equivalence margin supplied for omics: %s", paste(setdiff(omics, names(margin)), collapse = ", "))
  as.numeric(margin[omics])
}

.adjust_within_omic <- function(p, omic, method = "BH") {
  out <- rep(NA_real_, length(p))
  for (o in unique(omic)) {
    ii <- which(omic == o & is.finite(p))
    if (length(ii)) out[ii] <- stats::p.adjust(p[ii], method = method)
  }
  out
}

.pattern_from_values <- function(values, margins, trajectory_margin = 0.15) {
  if (length(values) < 2L || all(!is.finite(values))) return("insufficient")
  state <- ifelse(!is.finite(values), "uncertain",
                  ifelse(values > margins, "positive",
                         ifelse(values < -margins, "negative", "equivalent")))
  sigstates <- state[state %in% c("positive", "negative")]
  if (length(unique(sigstates)) > 1L) return("inversion")
  valid <- which(state != "uncertain")
  if (!length(valid)) return("insufficient")
  first_nonmissing <- valid[1]
  last_nonmissing <- utils::tail(valid, 1)
  if (state[first_nonmissing] == "equivalent" && first_nonmissing < length(state) &&
      any(state[seq.int(first_nonmissing + 1L, length(state))] %in% c("positive", "negative"), na.rm = TRUE)) return("emergence")
  if (state[first_nonmissing] %in% c("positive", "negative") && state[last_nonmissing] == "equivalent") return("buffering")
  if (length(sigstates) >= 2L && length(unique(sigstates)) == 1L) {
    ii <- which(state %in% c("positive", "negative"))
    if (length(ii) >= 2L) {
      aligned <- if (sigstates[1] == "negative") -values[ii] else values[ii]
      slope <- stats::coef(stats::lm(aligned ~ ii))[2]
      if (is.finite(slope) && slope <= -abs(trajectory_margin)) return("attenuation")
      if (is.finite(slope) && slope >=  abs(trajectory_margin)) return("amplification")
      return(if (sigstates[1] == "positive") "concordant_increase" else "concordant_decrease")
    }
  }
  if (all(state == "equivalent", na.rm = TRUE)) return("null_equivalent")
  "uncertain"
}

.hedges_g_se_matrix <- function(mat, group, reference, comparison, min_n = 3L) {
  i0 <- which(group == reference)
  i1 <- which(group == comparison)
  if (!length(i0) || !length(i1)) {
    n <- nrow(mat)
    return(list(g = rep(NA_real_, n), se = rep(NA_real_, n), z = rep(NA_real_, n)))
  }
  x0 <- mat[, i0, drop = FALSE]
  x1 <- mat[, i1, drop = FALSE]
  n0 <- rowSums(!is.na(x0)); n1 <- rowSums(!is.na(x1))
  m0 <- .row_mean_na(x0); m1 <- .row_mean_na(x1)
  v0 <- .row_var_na(x0); v1 <- .row_var_na(x1)
  df <- n0 + n1 - 2
  sp2 <- ((n0 - 1) * v0 + (n1 - 1) * v1) / df
  d <- (m1 - m0) / sqrt(sp2)
  g <- .hedges_J(df) * d
  varg <- (n0 + n1) / (n0 * n1) + (g * g) / (2 * pmax(n0 + n1, 1))
  se <- sqrt(varg)
  z <- g / se
  bad <- n0 < min_n | n1 < min_n | !is.finite(g) | !is.finite(se) | se <= 0 | sp2 <= 0
  g[bad] <- se[bad] <- z[bad] <- NA_real_
  list(g = g, se = se, z = z)
}
