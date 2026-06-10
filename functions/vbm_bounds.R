# ============================================================
# Core variance-based sensitivity utilities
# ============================================================

finite_sample_variance <- function(x) {

  x <- x[is.finite(x)]

  if (length(x) < 5) {
    return(NA_real_)
  }

  stats::var(x)
}

bias_scale_base <- function(w, Y, Z, config = NULL) {

  w0 <- w[Z == 0]
  Y0 <- Y[Z == 0]

  if (length(w0) < 5 || length(Y0) < 5) {
    return(NA_real_)
  }

  y_var <- finite_sample_variance(Y0)
  w_var <- finite_sample_variance(w0)

  if (!is.finite(y_var) || !is.finite(w_var) || y_var <= 0 || w_var <= 0) {
    return(NA_real_)
  }

  sqrt(y_var * w_var)
}

bias_correlation_limit <- function(w, Y, Z, config = NULL) {

  w0 <- w[Z == 0]
  Y0 <- Y[Z == 0]

  if (length(w0) < 5 || length(Y0) < 5) {
    return(NA_real_)
  }

  rho <- suppressWarnings(
    stats::cor(w0, Y0)
  )

  if (!is.finite(rho)) {
    return(NA_real_)
  }

  rho <- max(-1, min(1, rho))

  sqrt(max(0, 1 - rho^2))
}

bias_scale <- function(w, Y, Z, config = NULL) {

  base <- bias_scale_base(
    w = w,
    Y = Y,
    Z = Z,
    config = config
  )

  corr_limit <- bias_correlation_limit(
    w = w,
    Y = Y,
    Z = Z,
    config = config
  )

  if (!is.finite(base) || !is.finite(corr_limit)) {
    return(NA_real_)
  }

  base * corr_limit
}

bias_bound <- function(w, Y, Z, R2, config = NULL) {

  R2 <- max(0, min(R2, 1 - 1e-12))

  if (R2 == 0) {
    return(0)
  }

  scale <- bias_scale(
    w = w,
    Y = Y,
    Z = Z,
    config = config
  )

  if (!is.finite(scale)) {
    return(NA_real_)
  }

  scale * sqrt(R2 / (1 - R2))
}
