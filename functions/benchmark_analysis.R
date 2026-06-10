# ============================================================
# Covariate benchmark analysis
# ============================================================
#
# Base VBM and MSM benchmark rows are assembled here. The two Figure 3
# extension rows, "VBM, w/ Corr." and "MSM (Qbal)", are implemented in
# their own files and called from the shared benchmark workflow below.
# ============================================================

.config_value <- function(config, name, default) {
  if (is.null(config) || is.null(config[[name]])) {
    return(default)
  }
  config[[name]]
}

resolve_benchmark_groups <- function(config) {

  if (!is.null(config$benchmark_groups) && length(config$benchmark_groups) > 0) {
    return(config$benchmark_groups)
  }

  vars <- config$benchmark_vars
  if (is.null(vars) || length(vars) == 0) {
    stop("No benchmark variables or benchmark groups were specified.")
  }

  out <- as.list(vars)
  names(out) <- vars
  out
}

format_removed_terms <- function(remove_terms) {
  paste(remove_terms, collapse = " + ")
}

make_reduced_formula <- function(full_formula, remove_terms) {

  remove_terms <- remove_terms[!is.na(remove_terms) & nzchar(remove_terms)]

  if (length(remove_terms) == 0) {
    return(full_formula)
  }

  update(
    full_formula,
    paste(". ~ . -", paste(remove_terms, collapse = " - "))
  )
}

benchmark_group_index <- function(data, group) {

  if (is.null(group) || is.na(group)) {
    group <- "treated_weight_function"
  }

  if (group %in% c("treated", "treated_weight_function", "A1")) {
    return(data$Z == 1)
  }

  if (group %in% c("control", "control_weight_function", "Z0")) {
    return(data$Z == 0)
  }

  if (group %in% c("all", "all_weight_function")) {
    return(rep(TRUE, nrow(data)))
  }

  stop(paste0("Unknown benchmark group: ", group))
}

benchmark_sample_variance <- function(x) {

  x <- x[is.finite(x)]

  if (length(x) < 5) {
    return(NA_real_)
  }

  stats::var(x)
}

convert_r2_minus <- function(r2_minus_raw) {

  if (is.finite(r2_minus_raw) && r2_minus_raw > 0) {
    return(r2_minus_raw / (1 + r2_minus_raw))
  }

  0
}

compute_benchmark_R2_value <- function(
    w_full,
    w_reduced,
    data,
    group,
    config
) {

  idx <- benchmark_group_index(data, group)

  w_full_g <- w_full[idx]
  w_reduced_g <- w_reduced[idx]

  v_full <- benchmark_sample_variance(w_full_g)
  v_reduced <- benchmark_sample_variance(w_reduced_g)

  min_var <- .config_value(config, "benchmark_min_var", 1e-10)

  r2_minus_raw <- NA_real_
  r2_directional <- NA_real_
  status <- "paper_directional_R2"

  if (!is.finite(v_full) || !is.finite(v_reduced) ||
      v_full <= min_var || v_reduced <= min_var) {
    status <- "invalid_weight_variance"
    r2_used <- NA_real_
  } else {
    r2_minus_raw <- 1 - v_reduced / v_full
    r2_directional <- convert_r2_minus(r2_minus_raw)
    r2_used <- max(0, min(r2_directional, 1 - 1e-12))

    if (r2_minus_raw <= 0) {
      status <- "nonpositive_directional_R2_clamped_to_zero"
    }
  }

  list(
    group = group,
    R2_used = r2_used,
    R2_directional = r2_directional,
    R2_minus_raw = r2_minus_raw,
    status = status,
    var_full = v_full,
    var_reduced = v_reduced,
    w_full_q95 = as.numeric(quantile(w_full_g, probs = 0.95, na.rm = TRUE, names = FALSE)),
    w_full_q99 = as.numeric(quantile(w_full_g, probs = 0.99, na.rm = TRUE, names = FALSE)),
    w_reduced_q95 = as.numeric(quantile(w_reduced_g, probs = 0.95, na.rm = TRUE, names = FALSE)),
    w_reduced_q99 = as.numeric(quantile(w_reduced_g, probs = 0.99, na.rm = TRUE, names = FALSE)),
    cor_full_reduced_weights = suppressWarnings(cor(w_full_g, w_reduced_g))
  )
}

compute_benchmark_R2 <- function(
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

  selected_group <- .config_value(config, "benchmark_group", "treated_weight_function")

  selected <- compute_benchmark_R2_value(
    w_full = w_full,
    w_reduced = w_reduced,
    data = data,
    group = selected_group,
    config = config
  )

  treated <- compute_benchmark_R2_value(
    w_full = w_full,
    w_reduced = w_reduced,
    data = data,
    group = "treated_weight_function",
    config = config
  )

  control <- compute_benchmark_R2_value(
    w_full = w_full,
    w_reduced = w_reduced,
    data = data,
    group = "control_weight_function",
    config = config
  )

  all_units <- compute_benchmark_R2_value(
    w_full = w_full,
    w_reduced = w_reduced,
    data = data,
    group = "all_weight_function",
    config = config
  )

  full_cols <- colnames(model.matrix(full_formula, data))
  reduced_cols <- colnames(model.matrix(reduced_formula, data))

  list(
    benchmark_name = benchmark_name,
    removed_terms = format_removed_terms(remove_terms),
    reduced_formula = paste(deparse(reduced_formula), collapse = " "),
    R2_used = selected$R2_used,
    benchmark_status = selected$status,
    benchmark_group_used = selected_group,
    R2_directional = selected$R2_directional,
    R2_minus_raw = selected$R2_minus_raw,
    var_full = selected$var_full,
    var_reduced = selected$var_reduced,
    R2_control = control$R2_used,
    R2_treated = treated$R2_used,
    R2_all = all_units$R2_used,
    R2_control_directional = control$R2_directional,
    R2_treated_directional = treated$R2_directional,
    R2_all_directional = all_units$R2_directional,
    ps_full_min = min(ps_full, na.rm = TRUE),
    ps_full_max = max(ps_full, na.rm = TRUE),
    ps_reduced_min = min(ps_reduced, na.rm = TRUE),
    ps_reduced_max = max(ps_reduced, na.rm = TRUE),
    w_full_q95 = selected$w_full_q95,
    w_full_q99 = selected$w_full_q99,
    w_reduced_q95 = selected$w_reduced_q95,
    w_reduced_q99 = selected$w_reduced_q99,
    cor_full_reduced_weights = selected$cor_full_reduced_weights,
    full_model_n_columns = length(full_cols),
    reduced_model_n_columns = length(reduced_cols),
    removed_model_columns = paste(setdiff(full_cols, reduced_cols), collapse = " | "),
    new_model_columns = paste(setdiff(reduced_cols, full_cols), collapse = " | ")
  )
}

make_benchmark_row <- function(
    variable,
    method,
    sensitivity_parameter,
    sensitivity_value,
    lower,
    upper,
    removed_terms = NA_character_,
    interval_type = NA_character_,
    inference_type = NA_character_,
    tau_hat = NA_real_,
    deterministic_lower = NA_real_,
    deterministic_upper = NA_real_,
    bootstrap_lower = NA_real_,
    bootstrap_upper = NA_real_,
    bootstrap_B = NA_real_,
    bootstrap_n_valid = NA_real_,
    bootstrap_required_valid = NA_real_,
    bootstrap_pool = NA_character_,
    benchmark_status = NA_character_,
    benchmark_group_used = NA_character_,
    R2_minus_raw = NA_real_,
    R2_directional = NA_real_,
    R2_control = NA_real_,
    R2_treated = NA_real_,
    R2_all = NA_real_,
    R2_control_directional = NA_real_,
    R2_treated_directional = NA_real_,
    R2_all_directional = NA_real_,
    var_full = NA_real_,
    var_reduced = NA_real_,
    gamma_raw = NA_real_,
    gamma_used = NA_real_,
    gamma_group_used = NA_character_,
    ratio_max = NA_real_,
    ratio_q95 = NA_real_,
    ratio_q99 = NA_real_,
    ratio_qbal = NA_real_,
    qbal_effective_gamma = NA_real_,
    qbal_ratio_quantile = NA_real_,
    n_ratio = NA_real_,
    corr_raw = NA_real_,
    corr_used = NA_real_,
    corr_limit = NA_real_,
    corr_group_used = NA_character_,
    full_model_n_columns = NA_real_,
    reduced_model_n_columns = NA_real_,
    removed_model_columns = NA_character_
) {

  data.frame(
    variable = variable,
    removed_terms = removed_terms,
    method = method,
    sensitivity_parameter = sensitivity_parameter,
    sensitivity_value = sensitivity_value,
    lower = lower,
    upper = upper,
    interval_type = interval_type,
    inference_type = inference_type,
    tau_hat = tau_hat,
    deterministic_lower = deterministic_lower,
    deterministic_upper = deterministic_upper,
    bootstrap_lower = bootstrap_lower,
    bootstrap_upper = bootstrap_upper,
    bootstrap_B = bootstrap_B,
    bootstrap_n_valid = bootstrap_n_valid,
    bootstrap_required_valid = bootstrap_required_valid,
    bootstrap_pool = bootstrap_pool,
    benchmark_status = benchmark_status,
    benchmark_group_used = benchmark_group_used,
    R2_minus_raw = R2_minus_raw,
    R2_directional = R2_directional,
    R2_control = R2_control,
    R2_treated = R2_treated,
    R2_all = R2_all,
    R2_control_directional = R2_control_directional,
    R2_treated_directional = R2_treated_directional,
    R2_all_directional = R2_all_directional,
    var_full = var_full,
    var_reduced = var_reduced,
    gamma_raw = gamma_raw,
    gamma_used = gamma_used,
    gamma_group_used = gamma_group_used,
    ratio_max = ratio_max,
    ratio_q95 = ratio_q95,
    ratio_q99 = ratio_q99,
    ratio_qbal = ratio_qbal,
    qbal_effective_gamma = qbal_effective_gamma,
    qbal_ratio_quantile = qbal_ratio_quantile,
    n_ratio = n_ratio,
    corr_raw = corr_raw,
    corr_used = corr_used,
    corr_limit = corr_limit,
    corr_group_used = corr_group_used,
    full_model_n_columns = full_model_n_columns,
    reduced_model_n_columns = reduced_model_n_columns,
    removed_model_columns = removed_model_columns,
    stringsAsFactors = FALSE
  )
}

select_nearest_R2_row <- function(
    vbm_bootstrap_curve,
    R2,
    tolerance = 1e-10
) {

  if (is.null(vbm_bootstrap_curve) || nrow(vbm_bootstrap_curve) == 0) {
    return(NULL)
  }

  if (!is.finite(R2)) {
    return(NULL)
  }

  candidates <- vbm_bootstrap_curve[
    is.finite(vbm_bootstrap_curve$R2),
    ,
    drop = FALSE
  ]

  if (nrow(candidates) == 0) {
    return(NULL)
  }

  idx <- which.min(abs(candidates$R2 - R2))

  if (length(idx) == 0 || is.na(idx)) {
    return(NULL)
  }

  if (abs(candidates$R2[idx] - R2) > tolerance) {
    return(NULL)
  }

  candidates[idx, , drop = FALSE]
}

generate_vbm_benchmark_bootstrap_curve <- function(
    analysis,
    R2_values,
    config
) {

  enabled <- isTRUE(
    .config_value(config, "vbm_benchmark_bootstrap_enabled", TRUE)
  )

  if (!enabled) {
    return(data.frame())
  }

  R2_values <- sort(unique(R2_values[is.finite(R2_values)]))
  R2_values <- pmax(0, pmin(R2_values, 1 - 1e-12))

  if (length(R2_values) == 0) {
    return(data.frame())
  }

  B <- .config_value(
    config,
    "vbm_benchmark_bootstrap_B",
    .config_value(config, "fine_B", 300)
  )

  if (is.null(B) || !is.finite(B) || B <= 0) {
    return(data.frame())
  }

  seed_offset <- .config_value(
    config,
    "vbm_benchmark_bootstrap_seed_offset",
    8100
  )

  if (!is.null(config$seed) && is.finite(config$seed)) {
    set.seed(as.integer(config$seed + seed_offset))
  }

  stats <- generate_bootstrap_stats(
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

  for (r in R2_values) {

    ci <- bootstrap_ci_from_stats(
      stats = stats,
      R2 = r,
      B = B_requested,
      config = config
    )

    curve <- rbind(
      curve,
      data.frame(
        R2 = r,
        lower = ci["lower"],
        upper = ci["upper"],
        n_valid = ci["n_valid"],
        B_requested = B_requested,
        required_valid = required_valid,
        bootstrap_pool = "shared_fixed_benchmark_R2_pool",
        benchmark_bootstrap_type = "fixed_original_sample_R2",
        stringsAsFactors = FALSE
      )
    )
  }

  attr(curve, "bootstrap_stats") <- stats
  attr(curve, "B_requested") <- B_requested

  attr(curve, "bootstrap_stats_summary") <- bootstrap_stats_summary(
    stats = stats,
    stage = "benchmark_vbm",
    scope = "shared_fixed_benchmark_R2_pool"
  )

  curve
}

select_nearest_gamma_row <- function(
    msm_bootstrap_curve,
    gamma,
    tolerance = 1e-10
) {

  if (is.null(msm_bootstrap_curve) || nrow(msm_bootstrap_curve) == 0 ||
      !is.finite(gamma)) {
    return(NULL)
  }

  diffs <- abs(msm_bootstrap_curve$Gamma - gamma)
  if (all(!is.finite(diffs))) {
    return(NULL)
  }

  idx <- which.min(diffs)

  if (!is.finite(diffs[idx]) || diffs[idx] > tolerance) {
    warning(
      paste0(
        "No exact MSM bootstrap match for Gamma = ", gamma,
        "; using nearest Gamma = ", msm_bootstrap_curve$Gamma[idx], "."
      )
    )
  }

  msm_bootstrap_curve[idx, , drop = FALSE]
}

generate_msm_benchmark_bootstrap_curve <- function(
    analysis,
    gamma_values,
    config
) {

  enabled <- isTRUE(
    .config_value(config, "msm_benchmark_bootstrap_enabled", TRUE)
  )

  if (!enabled) {
    return(data.frame())
  }

  gamma_values <- sort(unique(gamma_values[is.finite(gamma_values) & gamma_values >= 1]))

  if (length(gamma_values) == 0) {
    return(data.frame())
  }

  B <- .config_value(
    config,
    "msm_benchmark_bootstrap_B",
    .config_value(config, "msm_bootstrap_B", 0)
  )

  if (is.null(B) || !is.finite(B) || B <= 0) {
    return(data.frame())
  }

  seed_offset <- .config_value(
    config,
    "msm_benchmark_bootstrap_seed_offset",
    9100
  )

  if (!is.null(config$seed) && is.finite(config$seed)) {
    set.seed(as.integer(config$seed + seed_offset))
  }

  curve <- generate_msm_bootstrap_curve(
    data = analysis,
    config = config,
    B = as.integer(B),
    Gamma_grid = gamma_values
  )

  if (nrow(curve) > 0) {
    curve$bootstrap_pool <- "shared_fixed_benchmark_Gamma_pool"
    curve$benchmark_bootstrap_type <- "fixed_original_sample_Gamma"
  }

  curve
}

generate_benchmark_comparison_data <- function(
    analysis,
    tau_hat,
    config,
    full_formula = ps_formula
) {

  out <- data.frame()
  benchmark_groups <- resolve_benchmark_groups(config)

  benchmark_items <- list()
  gamma_values <- c()
  R2_values <- c()

  include_qbal <- isTRUE(.config_value(config, "figure3_include_qbal", TRUE))
  include_vbm_corr <- isTRUE(.config_value(config, "figure3_include_vbm_corr", TRUE))

  for (v in names(benchmark_groups)) {

    remove_terms <- benchmark_groups[[v]]

    r2_info <- compute_benchmark_R2(
      full_formula = full_formula,
      benchmark_name = v,
      remove_terms = remove_terms,
      data = analysis,
      config = config
    )

    r2 <- r2_info$R2_used
    R2_values <- c(R2_values, r2)

    if (is.finite(r2)) {
      bb <- bias_bound(
        w = analysis$w,
        Y = analysis$Y,
        Z = analysis$Z,
        R2 = r2,
        config = config
      )
    } else {
      bb <- NA_real_
    }

    corr_info <- NULL
    bb_corr <- NA_real_

    if (include_vbm_corr) {
      corr_info <- compute_benchmark_corr(
        full_formula = full_formula,
        benchmark_name = v,
        remove_terms = remove_terms,
        data = analysis,
        config = config
      )

      if (is.finite(r2)) {
        bb_corr <- bias_bound_corr(
          w = analysis$w,
          Y = analysis$Y,
          Z = analysis$Z,
          R2 = r2,
          corr_bound = corr_info$corr_used,
          config = config
        )
      }
    }

    gamma_info <- compute_benchmark_Gamma(
      full_formula = full_formula,
      benchmark_name = v,
      remove_terms = remove_terms,
      data = analysis,
      config = config
    )

    gamma <- gamma_info$gamma_used
    gamma_values <- c(gamma_values, gamma)

    if (is.finite(gamma)) {
      msm_det_interval <- msm_att_bounds(
        analysis = analysis,
        Gamma = gamma,
        config = config
      )
    } else {
      msm_det_interval <- data.frame(
        lower = NA_real_,
        upper = NA_real_
      )
    }

    qbal_gamma <- NA_real_
    qbal_det_interval <- data.frame(
      lower = NA_real_,
      upper = NA_real_
    )

    if (include_qbal && is.finite(gamma)) {
      qbal_gamma <- qbal_effective_gamma_from_benchmark(
        gamma_info = gamma_info,
        config = config
      )
      gamma_values <- c(gamma_values, qbal_gamma)

      qbal_det_interval <- msm_qbal_att_bounds(
        analysis = analysis,
        Gamma = gamma,
        effective_Gamma = qbal_gamma,
        config = config
      )
    }

    benchmark_items[[v]] <- list(
      remove_terms = remove_terms,
      r2_info = r2_info,
      r2 = r2,
      bb = bb,
      corr_info = corr_info,
      bb_corr = bb_corr,
      gamma_info = gamma_info,
      gamma = gamma,
      msm_det_interval = msm_det_interval,
      qbal_gamma = qbal_gamma,
      qbal_det_interval = qbal_det_interval
    )
  }

  vbm_benchmark_bootstrap_curve <- generate_vbm_benchmark_bootstrap_curve(
    analysis = analysis,
    R2_values = R2_values,
    config = config
  )

  vbm_benchmark_bootstrap_stats <- attr(
    vbm_benchmark_bootstrap_curve,
    "bootstrap_stats"
  )

  vbm_benchmark_B_requested <- attr(
    vbm_benchmark_bootstrap_curve,
    "B_requested"
  )

  if (is.null(vbm_benchmark_B_requested) || !is.finite(vbm_benchmark_B_requested)) {
    vbm_benchmark_B_requested <- .config_value(
      config,
      "vbm_benchmark_bootstrap_B",
      .config_value(config, "fine_B", 300)
    )
  }

  msm_benchmark_bootstrap_curve <- generate_msm_benchmark_bootstrap_curve(
    analysis = analysis,
    gamma_values = gamma_values,
    config = config
  )

  for (v in names(benchmark_items)) {

    item <- benchmark_items[[v]]
    r2_info <- item$r2_info
    r2 <- item$r2
    bb <- item$bb
    corr_info <- item$corr_info
    bb_corr <- item$bb_corr
    gamma_info <- item$gamma_info
    gamma <- item$gamma
    msm_det_interval <- item$msm_det_interval
    qbal_gamma <- item$qbal_gamma
    qbal_det_interval <- item$qbal_det_interval
    remove_terms <- item$remove_terms

    vbm_boot_row <- select_nearest_R2_row(
      vbm_bootstrap_curve = vbm_benchmark_bootstrap_curve,
      R2 = r2
    )

    det_vbm_lower <- tau_hat - bb
    det_vbm_upper <- tau_hat + bb

    if (!is.null(vbm_boot_row) &&
        is.finite(vbm_boot_row$lower) &&
        is.finite(vbm_boot_row$upper)) {
      vbm_lower <- vbm_boot_row$lower
      vbm_upper <- vbm_boot_row$upper
      vbm_interval_type <- "bootstrap_percentile_fixed_R2"
      vbm_inference_type <- "bootstrap"
      vbm_bootstrap_lower <- vbm_boot_row$lower
      vbm_bootstrap_upper <- vbm_boot_row$upper
      vbm_bootstrap_B <- vbm_boot_row$B_requested
      vbm_bootstrap_n_valid <- vbm_boot_row$n_valid
      vbm_bootstrap_required_valid <- vbm_boot_row$required_valid
      vbm_bootstrap_pool <- vbm_boot_row$bootstrap_pool
      vbm_status <- paste(
        r2_info$benchmark_status,
        "bootstrap_fixed_original_R2",
        sep = " | "
      )
    } else {
      vbm_lower <- det_vbm_lower
      vbm_upper <- det_vbm_upper
      vbm_interval_type <- "closed_form_fallback_no_valid_benchmark_bootstrap"
      vbm_inference_type <- "closed_form_fallback"
      vbm_bootstrap_lower <- NA_real_
      vbm_bootstrap_upper <- NA_real_
      vbm_bootstrap_B <- NA_real_
      vbm_bootstrap_n_valid <- NA_real_
      vbm_bootstrap_required_valid <- NA_real_
      vbm_bootstrap_pool <- NA_character_
      vbm_status <- paste(
        r2_info$benchmark_status,
        "benchmark_bootstrap_unavailable_used_closed_form_bounds",
        sep = " | "
      )
    }

    out <- rbind(
      out,
      make_benchmark_row(
        variable = v,
        removed_terms = r2_info$removed_terms,
        method = "VBM benchmark (bootstrap)",
        sensitivity_parameter = "R2",
        sensitivity_value = r2,
        lower = vbm_lower,
        upper = vbm_upper,
        interval_type = vbm_interval_type,
        inference_type = vbm_inference_type,
        tau_hat = tau_hat,
        deterministic_lower = det_vbm_lower,
        deterministic_upper = det_vbm_upper,
        bootstrap_lower = vbm_bootstrap_lower,
        bootstrap_upper = vbm_bootstrap_upper,
        bootstrap_B = vbm_bootstrap_B,
        bootstrap_n_valid = vbm_bootstrap_n_valid,
        bootstrap_required_valid = vbm_bootstrap_required_valid,
        bootstrap_pool = vbm_bootstrap_pool,
        benchmark_status = vbm_status,
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
        full_model_n_columns = r2_info$full_model_n_columns,
        reduced_model_n_columns = r2_info$reduced_model_n_columns,
        removed_model_columns = r2_info$removed_model_columns
      )
    )

    if (include_vbm_corr) {
      out <- rbind(
        out,
        make_vbm_corr_benchmark_row(
          variable = v,
          tau_hat = tau_hat,
          r2_info = r2_info,
          r2 = r2,
          bb_corr = bb_corr,
          corr_info = corr_info,
          vbm_benchmark_bootstrap_stats = vbm_benchmark_bootstrap_stats,
          vbm_benchmark_B_requested = vbm_benchmark_B_requested,
          config = config
        )
      )
    }

    boot_row <- select_nearest_gamma_row(
      msm_bootstrap_curve = msm_benchmark_bootstrap_curve,
      gamma = gamma
    )

    if (!is.null(boot_row) &&
        is.finite(boot_row$lower) &&
        is.finite(boot_row$upper)) {
      msm_lower <- boot_row$lower
      msm_upper <- boot_row$upper
      interval_type <- "bootstrap_percentile_fixed_Gamma"
      inference_type <- "bootstrap"
      bootstrap_lower <- boot_row$lower
      bootstrap_upper <- boot_row$upper
      bootstrap_B <- boot_row$B_requested
      bootstrap_n_valid <- boot_row$n_valid
      bootstrap_required_valid <- boot_row$required_valid
      bootstrap_pool <- boot_row$bootstrap_pool
      msm_status <- paste(
        gamma_info$status,
        "bootstrap_fixed_original_Gamma",
        sep = " | "
      )
    } else {
      msm_lower <- msm_det_interval$lower
      msm_upper <- msm_det_interval$upper
      interval_type <- "closed_form_fallback_no_valid_benchmark_bootstrap"
      inference_type <- "closed_form_fallback"
      bootstrap_lower <- NA_real_
      bootstrap_upper <- NA_real_
      bootstrap_B <- NA_real_
      bootstrap_n_valid <- NA_real_
      bootstrap_required_valid <- NA_real_
      bootstrap_pool <- NA_character_
      msm_status <- paste(
        gamma_info$status,
        "benchmark_bootstrap_unavailable_used_closed_form_bounds",
        sep = " | "
      )
    }

    out <- rbind(
      out,
      make_benchmark_row(
        variable = v,
        removed_terms = format_removed_terms(remove_terms),
        method = "MSM benchmark (bootstrap)",
        sensitivity_parameter = "Gamma",
        sensitivity_value = gamma,
        lower = msm_lower,
        upper = msm_upper,
        interval_type = interval_type,
        inference_type = inference_type,
        tau_hat = tau_hat,
        deterministic_lower = msm_det_interval$lower,
        deterministic_upper = msm_det_interval$upper,
        bootstrap_lower = bootstrap_lower,
        bootstrap_upper = bootstrap_upper,
        bootstrap_B = bootstrap_B,
        bootstrap_n_valid = bootstrap_n_valid,
        bootstrap_required_valid = bootstrap_required_valid,
        bootstrap_pool = bootstrap_pool,
        benchmark_status = msm_status,
        gamma_raw = gamma_info$gamma_raw,
        gamma_used = gamma_info$gamma_used,
        gamma_group_used = gamma_info$gamma_group_used,
        ratio_max = gamma_info$ratio_max,
        ratio_q95 = gamma_info$ratio_q95,
        ratio_q99 = gamma_info$ratio_q99,
        ratio_qbal = gamma_info$ratio_qbal,
        n_ratio = gamma_info$n_ratio
      )
    )

    if (include_qbal) {
      out <- rbind(
        out,
        make_msm_qbal_benchmark_row(
          variable = v,
          remove_terms = remove_terms,
          tau_hat = tau_hat,
          gamma_info = gamma_info,
          gamma = gamma,
          qbal_gamma = qbal_gamma,
          qbal_det_interval = qbal_det_interval,
          msm_benchmark_bootstrap_curve = msm_benchmark_bootstrap_curve,
          config = config
        )
      )
    }
  }

  attr(out, "vbm_benchmark_bootstrap_curve") <- vbm_benchmark_bootstrap_curve
  attr(out, "msm_benchmark_bootstrap_curve") <- msm_benchmark_bootstrap_curve

  out
}

run_covariate_benchmarks <- function(
    analysis,
    tau_hat,
    full_formula = ps_formula,
    config
) {

  benchmark_df <- generate_benchmark_comparison_data(
    analysis = analysis,
    tau_hat = tau_hat,
    config = config,
    full_formula = full_formula
  )

  list(
    benchmark_df = benchmark_df,
    vbm_benchmark_bootstrap_curve = attr(benchmark_df, "vbm_benchmark_bootstrap_curve"),
    msm_benchmark_bootstrap_curve = attr(benchmark_df, "msm_benchmark_bootstrap_curve")
  )
}
