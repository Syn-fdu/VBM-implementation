# ============================================================
# VBM benchmark extension: outcome-correlation-calibrated bounds
# ============================================================
#
# This file contains only the Figure 3 benchmark variant
# "VBM, w/ Corr.".  The base VBM bound remains in
# vbm_bounds.R; this extension replaces the conservative
# correlation factor with an observed-covariate benchmark correlation.

bias_bound_corr <- function(w, Y, Z, R2, corr_bound, config = NULL) {

  R2 <- max(0, min(R2, 1 - 1e-12))

  if (R2 == 0) {
    return(0)
  }

  base <- bias_scale_base(
    w = w,
    Y = Y,
    Z = Z,
    config = config
  )

  if (!is.finite(base)) {
    return(NA_real_)
  }

  corr_limit <- bias_correlation_limit(
    w = w,
    Y = Y,
    Z = Z,
    config = config
  )

  if (!is.finite(corr_limit)) {
    corr_limit <- 1
  }

  corr_used <- abs(corr_bound)

  if (!is.finite(corr_used)) {
    corr_used <- corr_limit
  }

  corr_used <- max(0, min(corr_used, corr_limit, 1))

  base * corr_used * sqrt(R2 / (1 - R2))
}

compute_benchmark_corr <- function(
    full_formula,
    benchmark_name,
    remove_terms,
    data,
    config
) {

  bcfg <- make_benchmark_config(config)

  full_fit <- fit_ps_model(
    formula = full_formula,
    data = data,
    config = bcfg
  )

  ps_full <- clip_ps(
    ps = full_fit$fitted.values,
    lower = bcfg$ps_clip_lower,
    upper = bcfg$ps_clip_upper
  )

  reduced_formula <- make_reduced_formula(
    full_formula = full_formula,
    remove_terms = remove_terms
  )

  reduced_fit <- fit_ps_model(
    formula = reduced_formula,
    data = data,
    config = bcfg
  )

  ps_reduced <- clip_ps(
    ps = reduced_fit$fitted.values,
    lower = bcfg$ps_clip_lower,
    upper = bcfg$ps_clip_upper
  )

  w_full <- compute_att_weight_function(
    ps = ps_full,
    Z = data$Z
  )

  w_reduced <- compute_att_weight_function(
    ps = ps_reduced,
    Z = data$Z
  )

  selected_group <- .config_value(
    config,
    "vbm_corr_group",
    .config_value(config, "benchmark_group", "control_weight_function")
  )

  idx <- benchmark_group_index(data, selected_group)

  imbalance <- w_full[idx] - w_reduced[idx]
  y <- data$Y[idx]

  ok <- is.finite(imbalance) & is.finite(y)

  corr_raw <- NA_real_
  status <- "benchmark_weight_delta_outcome_correlation"

  if (sum(ok) >= 5 && stats::sd(imbalance[ok]) > 0 && stats::sd(y[ok]) > 0) {
    corr_raw <- suppressWarnings(stats::cor(imbalance[ok], y[ok]))
  } else {
    status <- "invalid_benchmark_correlation_set_to_zero"
  }

  corr_limit <- bias_correlation_limit(
    w = data$w,
    Y = data$Y,
    Z = data$Z,
    config = config
  )

  if (!is.finite(corr_limit)) {
    corr_limit <- 1
  }

  if (!is.finite(corr_raw)) {
    corr_used <- 0
  } else {
    corr_used <- min(abs(corr_raw), corr_limit)
  }

  list(
    benchmark_name = benchmark_name,
    removed_terms = format_removed_terms(remove_terms),
    corr_raw = corr_raw,
    corr_used = corr_used,
    corr_limit = corr_limit,
    corr_group_used = selected_group,
    status = status
  )
}

bootstrap_ci_from_stats_corr <- function(
    stats,
    R2,
    corr_bound,
    B,
    config
) {

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
        n_valid = n_valid
      )
    )
  }

  R2 <- max(0, min(R2, 1 - 1e-12))

  multiplier <- 0
  if (R2 > 0) {
    multiplier <- sqrt(R2 / (1 - R2))
  }

  corr_used <- abs(corr_bound)
  if (!is.finite(corr_used)) {
    corr_used <- 1
  }
  corr_used <- max(0, min(corr_used, 1))

  scale_vals <- stats$base_scale * corr_used

  lower_vals <- stats$tau - scale_vals * multiplier
  upper_vals <- stats$tau + scale_vals * multiplier

  c(
    lower = as.numeric(
      quantile(
        lower_vals,
        probs = config$percentile[1],
        na.rm = TRUE,
        names = FALSE
      )
    ),
    upper = as.numeric(
      quantile(
        upper_vals,
        probs = config$percentile[2],
        na.rm = TRUE,
        names = FALSE
      )
    ),
    n_valid = n_valid
  )
}

make_vbm_corr_benchmark_row <- function(
    variable,
    tau_hat,
    r2_info,
    r2,
    bb_corr,
    corr_info,
    vbm_benchmark_bootstrap_stats,
    vbm_benchmark_B_requested,
    config
) {

  det_lower <- tau_hat - bb_corr
  det_upper <- tau_hat + bb_corr

  ci <- c(lower = NA_real_, upper = NA_real_, n_valid = NA_real_)

  if (!is.null(vbm_benchmark_bootstrap_stats) &&
      nrow(vbm_benchmark_bootstrap_stats) > 0 &&
      !is.null(corr_info) && is.finite(r2)) {
    ci <- bootstrap_ci_from_stats_corr(
      stats = vbm_benchmark_bootstrap_stats,
      R2 = r2,
      corr_bound = corr_info$corr_used,
      B = vbm_benchmark_B_requested,
      config = config
    )
  }

  if (is.finite(ci["lower"]) && is.finite(ci["upper"])) {
    lower <- ci["lower"]
    upper <- ci["upper"]
    interval_type <- "bootstrap_percentile_fixed_R2_corr"
    inference_type <- "bootstrap"
    bootstrap_lower <- ci["lower"]
    bootstrap_upper <- ci["upper"]
    bootstrap_B <- vbm_benchmark_B_requested
    bootstrap_n_valid <- ci["n_valid"]
    bootstrap_required <- bootstrap_required_valid(
      B = vbm_benchmark_B_requested,
      config = config
    )
    bootstrap_pool <- "shared_fixed_benchmark_R2_corr_pool"
    status <- paste(
      corr_info$status,
      "bootstrap_fixed_original_R2_and_corr",
      sep = " | "
    )
  } else {
    lower <- det_lower
    upper <- det_upper
    interval_type <- "closed_form_corr_fallback_no_valid_benchmark_bootstrap"
    inference_type <- "closed_form_fallback"
    bootstrap_lower <- NA_real_
    bootstrap_upper <- NA_real_
    bootstrap_B <- NA_real_
    bootstrap_n_valid <- NA_real_
    bootstrap_required <- NA_real_
    bootstrap_pool <- NA_character_
    status <- paste(
      corr_info$status,
      "benchmark_bootstrap_unavailable_used_closed_form_corr_bounds",
      sep = " | "
    )
  }

  make_benchmark_row(
    variable = variable,
    removed_terms = r2_info$removed_terms,
    method = "VBM, w/ Corr. benchmark (bootstrap)",
    sensitivity_parameter = "R2 + Corr",
    sensitivity_value = r2,
    lower = lower,
    upper = upper,
    interval_type = interval_type,
    inference_type = inference_type,
    tau_hat = tau_hat,
    deterministic_lower = det_lower,
    deterministic_upper = det_upper,
    bootstrap_lower = bootstrap_lower,
    bootstrap_upper = bootstrap_upper,
    bootstrap_B = bootstrap_B,
    bootstrap_n_valid = bootstrap_n_valid,
    bootstrap_required_valid = bootstrap_required,
    bootstrap_pool = bootstrap_pool,
    benchmark_status = status,
    benchmark_group_used = r2_info$benchmark_group_used,
    R2_minus_raw = r2_info$R2_minus_raw,
    R2_directional = r2_info$R2_directional,
    R2_control = r2_info$R2_control,
    R2_treated = r2_info$R2_treated,
    R2_all = r2_info$R2_all,
    R2_control_directional = r2_info$R2_control_directional,
    R2_treated_directional = r2_info$R2_treated_directional,
    R2_all_directional = r2_info$R2_all_directional,
    var_full = r2_info$var_full,
    var_reduced = r2_info$var_reduced,
    corr_raw = corr_info$corr_raw,
    corr_used = corr_info$corr_used,
    corr_limit = corr_info$corr_limit,
    corr_group_used = corr_info$corr_group_used,
    full_model_n_columns = r2_info$full_model_n_columns,
    reduced_model_n_columns = r2_info$reduced_model_n_columns,
    removed_model_columns = r2_info$removed_model_columns
  )
}
