# ============================================================
# RGM conservative range-bound calculations for the integrated main.R workflow
# ============================================================
#
# The conservative RGM model uses the total-variation/L1 sensitivity
# parameter
#
#   T = 1/2 * sum_i |q_i - p_i|,
#
# where p is the observed normalized control-weight distribution and q is an
# admissible ideal control-weight distribution.  If the outcome is bounded by
# [Y_L, Y_U], then the ATT bias is bounded by
#
#   Bias_RGM(T) = T * (Y_U - Y_L).
#
# This file deliberately contains only the conservative closed-form/range-bound
# calculations.  It does not search for T*.  The sharp finite-sample RGM model
# is implemented separately in rgm_model/rgm_sharp.R.
# ============================================================

rgm_ext_value <- function(config, name, default) {
  if (is.null(config) || is.null(config[[name]])) {
    return(default)
  }
  config[[name]]
}

rgm_clamp_T <- function(T) {
  pmax(0, pmin(T, 1))
}

rgm_outcome_range <- function(Y, Z = NULL, config = NULL) {

  mode <- rgm_ext_value(config, "rgm_conservative_outcome_bound_mode", NA_character_)
  if (is.na(mode) || !nzchar(mode)) {
    mode <- rgm_ext_value(config, "rgm_outcome_bound_mode", "sample_control")
  }

  probs <- rgm_ext_value(config, "rgm_winsor_probs", c(0.01, 0.99))

  if (mode %in% c("known", "fixed")) {
    y_lower <- rgm_ext_value(config, "rgm_y_lower", NA_real_)
    y_upper <- rgm_ext_value(config, "rgm_y_upper", NA_real_)

    if (!is.finite(y_lower) || !is.finite(y_upper) || y_lower >= y_upper) {
      stop("For known/fixed conservative RGM bounds, set finite rgm_y_lower < rgm_y_upper.")
    }

    return(
      list(
        lower = y_lower,
        upper = y_upper,
        range = y_upper - y_lower,
        mode = mode,
        n_used = NA_integer_
      )
    )
  }

  y <- Y
  if (!is.null(Z) && mode %in% c("sample_control", "winsorized_control", "control")) {
    y <- Y[Z == 0]
  }

  y <- y[is.finite(y)]

  if (length(y) == 0) {
    return(
      list(
        lower = NA_real_,
        upper = NA_real_,
        range = NA_real_,
        mode = mode,
        n_used = 0L
      )
    )
  }

  if (mode %in% c("winsorized_control", "winsorized_all", "winsorized")) {
    probs <- as.numeric(probs)
    if (length(probs) != 2 || any(!is.finite(probs)) ||
        probs[1] < 0 || probs[2] > 1 || probs[1] >= probs[2]) {
      stop("rgm_winsor_probs must be c(lower_prob, upper_prob) inside [0, 1].")
    }

    bounds <- stats::quantile(
      y,
      probs = probs,
      na.rm = TRUE,
      names = FALSE,
      type = 7
    )
    y_lower <- as.numeric(bounds[1])
    y_upper <- as.numeric(bounds[2])
  } else {
    y_lower <- min(y, na.rm = TRUE)
    y_upper <- max(y, na.rm = TRUE)
  }

  list(
    lower = y_lower,
    upper = y_upper,
    range = y_upper - y_lower,
    mode = mode,
    n_used = length(y)
  )
}

rgm_conservative_bias_bound <- function(Y, Z, T, config = NULL) {

  T <- rgm_clamp_T(T)

  yr <- rgm_outcome_range(
    Y = Y,
    Z = Z,
    config = config
  )

  if (!is.finite(yr$range) || yr$range < 0) {
    return(NA_real_)
  }

  T * yr$range
}

rgm_conservative_att_bounds <- function(analysis, T, config = NULL) {

  if (!all(c("Y", "Z", "w") %in% names(analysis))) {
    stop("analysis must contain Y, Z, and w columns.")
  }

  T <- rgm_clamp_T(T)

  tau_hat <- tryCatch(
    estimate_att(
      Y = analysis$Y,
      Z = analysis$Z,
      w = analysis$w
    ),
    error = function(e) NA_real_
  )

  yr <- rgm_outcome_range(
    Y = analysis$Y,
    Z = analysis$Z,
    config = config
  )

  bias <- rgm_conservative_bias_bound(
    Y = analysis$Y,
    Z = analysis$Z,
    T = T,
    config = config
  )

  if (!is.finite(tau_hat) || !is.finite(bias)) {
    lower <- NA_real_
    upper <- NA_real_
    status <- "invalid_tau_or_outcome_range"
  } else {
    lower <- tau_hat - bias
    upper <- tau_hat + bias
    status <- "conservative_range_bound_closed_form"
  }

  data.frame(
    method = "RGM-conservative",
    bound_type = "conservative_range",
    T = T,
    lower = lower,
    upper = upper,
    tau_hat = tau_hat,
    bias = bias,
    outcome_lower = yr$lower,
    outcome_upper = yr$upper,
    outcome_range = yr$range,
    outcome_bound_mode = yr$mode,
    n_range_used = yr$n_used,
    status = status,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# Conservative benchmark bootstrap at fixed original-sample T values.
# This is not a T* search; it is only used for benchmark comparison rows.
# ------------------------------------------------------------

generate_rgm_conservative_bootstrap_stats <- function(data, B = 100, config) {

  tau_vals <- rep(NA_real_, B)
  range_vals <- rep(NA_real_, B)
  y_lower_vals <- rep(NA_real_, B)
  y_upper_vals <- rep(NA_real_, B)

  failures <- c(
    one_class = 0,
    ps_fit = 0,
    weights = 0,
    tau_or_range = 0
  )

  for (b in seq_len(B)) {

    d <- bootstrap_resample_data(
      data = data,
      config = config
    )

    if (length(unique(d$Z)) < 2) {
      failures["one_class"] <- failures["one_class"] + 1
      next
    }

    fit_b <- tryCatch(
      fit_ps_model(
        formula = ps_formula,
        data = d,
        config = config
      ),
      error = function(e) NULL
    )

    if (is.null(fit_b)) {
      failures["ps_fit"] <- failures["ps_fit"] + 1
      next
    }

    ps_b <- clip_ps(
      ps = fit_b$fitted.values,
      lower = config$ps_clip_lower,
      upper = config$ps_clip_upper
    )

    w_b <- compute_weights(
      ps = ps_b,
      Z = d$Z,
      config = config
    )

    if (any(!is.finite(w_b))) {
      failures["weights"] <- failures["weights"] + 1
      next
    }

    tau_b <- tryCatch(
      estimate_att(
        Y = d$Y,
        Z = d$Z,
        w = w_b
      ),
      error = function(e) NA_real_
    )

    yr_b <- rgm_outcome_range(
      Y = d$Y,
      Z = d$Z,
      config = config
    )

    if (!is.finite(tau_b) || !is.finite(yr_b$range)) {
      failures["tau_or_range"] <- failures["tau_or_range"] + 1
      next
    }

    tau_vals[b] <- tau_b
    range_vals[b] <- yr_b$range
    y_lower_vals[b] <- yr_b$lower
    y_upper_vals[b] <- yr_b$upper
  }

  stats <- data.frame(
    b = seq_len(B),
    tau = tau_vals,
    outcome_range = range_vals,
    outcome_lower = y_lower_vals,
    outcome_upper = y_upper_vals,
    stringsAsFactors = FALSE
  )

  stats <- stats[
    is.finite(stats$tau) & is.finite(stats$outcome_range),
    ,
    drop = FALSE
  ]

  attr(stats, "B_requested") <- B
  attr(stats, "failures") <- failures

  stats
}

rgm_conservative_bootstrap_ci_from_stats <- function(stats, T, B, config) {

  n_valid <- nrow(stats)
  required_valid <- bootstrap_required_valid(
    B = B,
    config = config
  )

  if (n_valid <= 0 || n_valid < required_valid) {
    return(
      c(
        lower = NA_real_,
        upper = NA_real_,
        n_valid = n_valid,
        required_valid = required_valid
      )
    )
  }

  T <- rgm_clamp_T(T)

  lower_vals <- stats$tau - T * stats$outcome_range
  upper_vals <- stats$tau + T * stats$outcome_range

  c(
    lower = as.numeric(
      stats::quantile(
        lower_vals,
        probs = config$percentile[1],
        na.rm = TRUE,
        names = FALSE
      )
    ),
    upper = as.numeric(
      stats::quantile(
        upper_vals,
        probs = config$percentile[2],
        na.rm = TRUE,
        names = FALSE
      )
    ),
    n_valid = n_valid,
    required_valid = required_valid
  )
}

generate_rgm_conservative_benchmark_bootstrap_curve <- function(
    analysis,
    T_values,
    config
) {

  enabled <- isTRUE(
    rgm_ext_value(config, "rgm_benchmark_bootstrap_enabled", TRUE)
  )

  if (!enabled) {
    return(data.frame())
  }

  T_values <- sort(unique(T_values[is.finite(T_values)]))
  T_values <- rgm_clamp_T(T_values)

  if (length(T_values) == 0) {
    return(data.frame())
  }

  B <- rgm_ext_value(
    config,
    "rgm_benchmark_bootstrap_B",
    rgm_ext_value(config, "fine_B", 300)
  )

  if (is.null(B) || !is.finite(B) || B <= 0) {
    return(data.frame())
  }

  seed_offset <- rgm_ext_value(
    config,
    "rgm_conservative_benchmark_seed_offset",
    11100
  )

  if (!is.null(config$seed) && is.finite(config$seed)) {
    set.seed(as.integer(config$seed + seed_offset))
  }

  stats <- generate_rgm_conservative_bootstrap_stats(
    data = analysis,
    B = as.integer(B),
    config = config
  )

  B_requested <- as.integer(B)
  required_valid <- bootstrap_required_valid(
    B = B_requested,
    config = config
  )

  curve <- data.frame()

  for (t in T_values) {

    ci <- rgm_conservative_bootstrap_ci_from_stats(
      stats = stats,
      T = t,
      B = B_requested,
      config = config
    )

    curve <- rbind(
      curve,
      data.frame(
        method = "RGM-conservative",
        bound_type = "conservative_range",
        T = t,
        lower = ci["lower"],
        upper = ci["upper"],
        n_valid = ci["n_valid"],
        B_requested = B_requested,
        required_valid = required_valid,
        bootstrap_pool = "shared_fixed_benchmark_T_pool_conservative",
        benchmark_bootstrap_type = "fixed_original_sample_T_conservative",
        stringsAsFactors = FALSE
      )
    )
  }

  attr(curve, "bootstrap_stats") <- stats
  attr(curve, "B_requested") <- B_requested
  attr(curve, "failures") <- attr(stats, "failures")

  curve
}
