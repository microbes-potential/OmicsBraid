#' Harmonize the sign orientation of omic-layer effects
#'
#' Some omic measurements have an interpretation whose natural direction is
#' opposite to an activity/abundance scale used in other layers. This helper
#' multiplies selected omic effects by +1 or -1 and applies the corresponding
#' sign transformation to covariance matrices. Use it only when the orientation
#' is scientifically justified; it must not be used to force apparent agreement.
#'
#' @param effects Effect table with `omic` and `effect` columns.
#' @param covariance Optional OmicsBraid covariance object or named covariance list.
#' @param orientation Named numeric vector with one value (+1 or -1) per omic to be re-oriented.
#'   Omics not named default to +1.
#' @return A list with oriented `effects` and `covariance`.
#' @export
orient_omics <- function(effects, covariance = NULL, orientation) {
  if (missing(orientation) || is.null(orientation)) return(list(effects = effects, covariance = covariance))
  if (!is.numeric(orientation) || is.null(names(orientation)) || any(!orientation %in% c(-1, 1))) {
    .stopf("'orientation' must be a named numeric vector containing only +1 and -1.")
  }
  if (anyDuplicated(names(orientation))) .stopf("'orientation' names must be unique omic names.")
  mult <- rep(1, nrow(effects))
  hit <- effects$omic %in% names(orientation)
  mult[hit] <- orientation[effects$omic[hit]]
  out <- effects
  old_effect <- out$effect
  out$effect <- old_effect * mult
  if (all(c("conf_low", "conf_high") %in% names(out))) {
    old_low <- out$conf_low; old_high <- out$conf_high
    flip <- mult < 0
    out$conf_low[flip] <- -old_high[flip]
    out$conf_high[flip] <- -old_low[flip]
  }
  attr(out, "orientation") <- orientation

  cv <- covariance
  transform_matrix <- function(M) {
    os <- rownames(M)
    if (is.null(os) || is.null(colnames(M))) .stopf("Covariance matrices must have omic row/column names for orientation.")
    m <- rep(1, length(os)); names(m) <- os
    hh <- os %in% names(orientation); m[hh] <- orientation[os[hh]]
    D <- diag(m, nrow = length(m), ncol = length(m))
    dimnames(D) <- list(os, os)
    D %*% M %*% D
  }
  if (inherits(cv, "omics_braid_covariance")) {
    cv$covariance <- lapply(cv$covariance, transform_matrix)
    cv$correlation <- lapply(cv$correlation, transform_matrix)
    if (!is.null(cv$boot_effects)) {
      for (o in intersect(names(cv$boot_effects), names(orientation))) cv$boot_effects[[o]] <- cv$boot_effects[[o]] * orientation[o]
    }
    cv$orientation <- orientation
  } else if (is.list(cv) && !is.null(cv)) {
    cv <- lapply(cv, transform_matrix)
  }
  list(effects = out, covariance = cv)
}
