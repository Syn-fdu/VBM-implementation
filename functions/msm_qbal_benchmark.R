# ============================================================
# MSM benchmark extension: quantile-ratio Qbal proxy
# ============================================================
#
# This file contains only the Figure 3 benchmark variant "MSM (Qbal)".
# The standard MSM bound remains in msm_analysis.R.  The Qbal row keeps
# the benchmark Gamma reported for MSM and computes the interval at an
# effective Gamma based on a high quantile of the full-vs-reduced log-ratio.

msm_qbal_ratio_summary <- function(log_abs_ratio, config) {

  log_abs_ratio <- log_abs_ratio[is.finite(log_abs_ratio)]

  if (length(log_abs_ratio) == 0) {
    return(
      list(
        ratio_q95 = NA_real_,
        ratio_q99 = NA_real_,
        ratio_qbal = NA_real_
      )
    )
  }

  qbal_q <- 0.95
  if (!is.null(config$msm_qbal_ratio_quantile) &&
      is.finite(config$msm_qbal_ratio_quantile)) {
    qbal_q <- config$msm_qbal_ratio_quantile
  }
  qbal_q <- max(0.50, min(qbal_q, 1.00))

  list(
    ratio_q95 = exp(stats::quantile(
      log_abs_ratio,
      probs = 0.95,
      na.rm = TRUE,
      names = FALSE
    )),
    ratio_q99 = exp(stats::quantile(
      log_abs_ratio,
      probs = 0.99,
      na.rm = TRUE,
      names = FALSE
    )),
    ratio_qbal = exp(stats::quantile(
      log_abs_ratio,
      probs = qbal_q,
      na.rm = TRUE,
      names = FALSE
    ))
  )
}

msm_qbal_att_bounds <- function(
    analysis,
    Gamma,
    effective_Gamma = NULL,
    config = NULL
) {

  if (is.null(effective_Gamma) || !is.finite(effective_Gamma)) {
    effective_Gamma <- Gamma
  }

  effective_Gamma <- max(1, min(effective_Gamma, Gamma))

  out <- msm_att_bounds(
    analysis = analysis,
    Gamma = effective_Gamma,
    config = config
  )

  out$Gamma <- Gamma
  out$qbal_effective_Gamma <- effective_Gamma
  out$status <- paste(
    "qbal_ratio_quantile_proxy",
    out$status,
    sep = " | "
  )

  out
}

qbal_effective_gamma_from_benchmark <- function(
    gamma_info,
    config
) {

  gamma <- gamma_info$gamma_used

  if (!is.finite(gamma) || gamma < 1) {
    return(NA_real_)
  }

  q <- 0.95
  if (!is.null(config$msm_qbal_ratio_quantile) &&
      is.finite(config$msm_qbal_ratio_quantile)) {
    q <- config$msm_qbal_ratio_quantile
  }

  q <- max(0.50, min(q, 1.00))

  if (is.finite(gamma_info$ratio_qbal)) {
    gamma_q <- gamma_info$ratio_qbal
  } else if (q <= 0.95 && is.finite(gamma_info$ratio_q95)) {
    gamma_q <- gamma_info$ratio_q95
  } else if (q <= 0.99 && is.finite(gamma_info$ratio_q99)) {
    gamma_q <- gamma_info$ratio_q99
  } else {
    gamma_q <- gamma
  }

  max(1, min(gamma, gamma_q))
}

make_msm_qbal_benchmark_row <- function(
    variable,
    remove_terms,
    tau_hat,
    gamma_info,
    gamma,
    qbal_gamma,
    qbal_det_interval,
    msm_benchmark_bootstrap_curve,
    config
) {

  qbal_boot_row <- select_nearest_gamma_row(
    msm_bootstrap_curve = msm_benchmark_bootstrap_curve,
    gamma = qbal_gamma
  )

  if (!is.null(qbal_boot_row) &&
      is.finite(qbal_boot_row$lower) &&
      is.finite(qbal_boot_row$upper)) {
    lower <- qbal_boot_row$lower
    upper <- qbal_boot_row$upper
    interval_type <- "bootstrap_percentile_fixed_qbal_Gamma"
    inference_type <- "bootstrap"
    bootstrap_lower <- qbal_boot_row$lower
    bootstrap_upper <- qbal_boot_row$upper
    bootstrap_B <- qbal_boot_row$B_requested
    bootstrap_n_valid <- qbal_boot_row$n_valid
    bootstrap_required <- qbal_boot_row$required_valid
    bootstrap_pool <- qbal_boot_row$bootstrap_pool
    status <- paste(
      gamma_info$status,
      "qbal_ratio_quantile_proxy",
      "bootstrap_fixed_original_effective_Gamma",
      sep = " | "
    )
  } else {
    lower <- qbal_det_interval$lower
    upper <- qbal_det_interval$upper
    interval_type <- "closed_form_qbal_fallback_no_valid_benchmark_bootstrap"
    inference_type <- "closed_form_fallback"
    bootstrap_lower <- NA_real_
    bootstrap_upper <- NA_real_
    bootstrap_B <- NA_real_
    bootstrap_n_valid <- NA_real_
    bootstrap_required <- NA_real_
    bootstrap_pool <- NA_character_
    status <- paste(
      gamma_info$status,
      "qbal_ratio_quantile_proxy",
      "benchmark_bootstrap_unavailable_used_closed_form_bounds",
      sep = " | "
    )
  }

  make_benchmark_row(
    variable = variable,
    removed_terms = format_removed_terms(remove_terms),
    method = "MSM (Qbal) benchmark (bootstrap)",
    sensitivity_parameter = "Gamma + Qbal",
    sensitivity_value = gamma,
    lower = lower,
    upper = upper,
    interval_type = interval_type,
    inference_type = inference_type,
    tau_hat = tau_hat,
    deterministic_lower = qbal_det_interval$lower,
    deterministic_upper = qbal_det_interval$upper,
    bootstrap_lower = bootstrap_lower,
    bootstrap_upper = bootstrap_upper,
    bootstrap_B = bootstrap_B,
    bootstrap_n_valid = bootstrap_n_valid,
    bootstrap_required_valid = bootstrap_required,
    bootstrap_pool = bootstrap_pool,
    benchmark_status = status,
    gamma_raw = gamma_info$gamma_raw,
    gamma_used = gamma_info$gamma_used,
    gamma_group_used = gamma_info$gamma_group_used,
    ratio_max = gamma_info$ratio_max,
    ratio_q95 = gamma_info$ratio_q95,
    ratio_q99 = gamma_info$ratio_q99,
    ratio_qbal = gamma_info$ratio_qbal,
    qbal_effective_gamma = qbal_gamma,
    qbal_ratio_quantile = .config_value(config, "msm_qbal_ratio_quantile", 0.95),
    n_ratio = gamma_info$n_ratio
  )
}
