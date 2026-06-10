# ============================================================
# RGM benchmark calculations and MSM/VBM/RGM comparison plot for the integrated main.R workflow
# ============================================================
#
# This file contains benchmark-specific RGM code used by the optional
# integrated extension block in main.R.  It reuses the original project's
# VBM/MSM benchmark functions and adds two RGM rows:
#   - RGM-sharp: finite-sample TV/L1 sharp bound
#   - RGM-conservative: closed-form range bound
# ============================================================

rgm_compute_benchmark_T_value <- function(
    w_full,
    w_reduced,
    data,
    group,
    config
) {

  idx <- benchmark_group_index(data, group)

  w_full_g <- w_full[idx]
  w_reduced_g <- w_reduced[idx]

  ok <- is.finite(w_full_g) & is.finite(w_reduced_g) & w_full_g > 0 & w_reduced_g > 0
  n_valid <- sum(ok)

  if (n_valid < 2) {
    return(
      list(
        group = group,
        T_used = NA_real_,
        l1_distance = NA_real_,
        status = "invalid_weight_probability_vectors",
        n_valid = n_valid,
        full_weight_mass = sum(w_full_g[is.finite(w_full_g) & w_full_g > 0]),
        reduced_weight_mass = sum(w_reduced_g[is.finite(w_reduced_g) & w_reduced_g > 0]),
        max_abs_probability_shift = NA_real_,
        q95_abs_probability_shift = NA_real_,
        cor_full_reduced_prob = NA_real_
      )
    )
  }

  w_full_g <- w_full_g[ok]
  w_reduced_g <- w_reduced_g[ok]

  p_full <- w_full_g / sum(w_full_g)
  p_reduced <- w_reduced_g / sum(w_reduced_g)

  l1 <- sum(abs(p_full - p_reduced))
  T <- rgm_clamp_T(0.5 * l1)

  list(
    group = group,
    T_used = T,
    l1_distance = l1,
    status = "normalized_full_vs_reduced_weight_TV_distance",
    n_valid = n_valid,
    full_weight_mass = sum(w_full_g),
    reduced_weight_mass = sum(w_reduced_g),
    max_abs_probability_shift = max(abs(p_full - p_reduced)),
    q95_abs_probability_shift = as.numeric(
      stats::quantile(
        abs(p_full - p_reduced),
        0.95,
        na.rm = TRUE,
        names = FALSE
      )
    ),
    cor_full_reduced_prob = suppressWarnings(stats::cor(p_full, p_reduced))
  )
}

rgm_compute_benchmark_T <- function(
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

  selected_group <- rgm_ext_value(
    config,
    "rgm_benchmark_group",
    rgm_ext_value(config, "benchmark_group", "control_weight_function")
  )

  selected <- rgm_compute_benchmark_T_value(
    w_full = w_full,
    w_reduced = w_reduced,
    data = data,
    group = selected_group,
    config = config
  )

  treated <- rgm_compute_benchmark_T_value(
    w_full = w_full,
    w_reduced = w_reduced,
    data = data,
    group = "treated_weight_function",
    config = config
  )

  control <- rgm_compute_benchmark_T_value(
    w_full = w_full,
    w_reduced = w_reduced,
    data = data,
    group = "control_weight_function",
    config = config
  )

  all_units <- rgm_compute_benchmark_T_value(
    w_full = w_full,
    w_reduced = w_reduced,
    data = data,
    group = "all_weight_function",
    config = config
  )

  full_cols <- colnames(stats::model.matrix(full_formula, data))
  reduced_cols <- colnames(stats::model.matrix(reduced_formula, data))

  list(
    benchmark_name = benchmark_name,
    removed_terms = format_removed_terms(remove_terms),
    reduced_formula = paste(deparse(reduced_formula), collapse = " "),
    T_used = selected$T_used,
    benchmark_status = selected$status,
    benchmark_group_used = selected_group,
    l1_distance = selected$l1_distance,
    T_control = control$T_used,
    T_treated = treated$T_used,
    T_all = all_units$T_used,
    l1_control = control$l1_distance,
    l1_treated = treated$l1_distance,
    l1_all = all_units$l1_distance,
    n_T = selected$n_valid,
    full_weight_mass = selected$full_weight_mass,
    reduced_weight_mass = selected$reduced_weight_mass,
    max_abs_probability_shift = selected$max_abs_probability_shift,
    q95_abs_probability_shift = selected$q95_abs_probability_shift,
    cor_full_reduced_prob = selected$cor_full_reduced_prob,
    full_model_n_columns = length(full_cols),
    reduced_model_n_columns = length(reduced_cols),
    removed_model_columns = paste(setdiff(full_cols, reduced_cols), collapse = " | "),
    new_model_columns = paste(setdiff(reduced_cols, full_cols), collapse = " | ")
  )
}

rgm_select_nearest_T_row <- function(curve, T, tolerance = 1e-10) {

  if (is.null(curve) || nrow(curve) == 0 || !is.finite(T)) {
    return(NULL)
  }

  diffs <- abs(curve$T - T)
  if (all(!is.finite(diffs))) {
    return(NULL)
  }

  idx <- which.min(diffs)

  if (!is.finite(diffs[idx]) || diffs[idx] > tolerance) {
    warning(
      paste0(
        "No exact RGM bootstrap match for T = ", T,
        "; using nearest T = ", curve$T[idx], "."
      )
    )
  }

  curve[idx, , drop = FALSE]
}

make_extension_benchmark_row <- function(
    variable,
    removed_terms,
    method,
    sensitivity_parameter,
    sensitivity_value,
    lower,
    upper,
    tau_hat,
    interval_type,
    inference_type,
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
    R2 = NA_real_,
    Gamma = NA_real_,
    T = NA_real_,
    R2_minus_raw = NA_real_,
    R2_directional = NA_real_,
    var_full = NA_real_,
    var_reduced = NA_real_,
    gamma_raw = NA_real_,
    ratio_max = NA_real_,
    ratio_q95 = NA_real_,
    ratio_q99 = NA_real_,
    n_ratio = NA_real_,
    l1_distance = NA_real_,
    T_control = NA_real_,
    T_treated = NA_real_,
    T_all = NA_real_,
    n_T = NA_real_,
    full_weight_mass_T = NA_real_,
    reduced_weight_mass_T = NA_real_,
    max_abs_probability_shift = NA_real_,
    q95_abs_probability_shift = NA_real_,
    cor_full_reduced_prob = NA_real_
) {

  data.frame(
    variable = variable,
    removed_terms = removed_terms,
    method = method,
    sensitivity_parameter = sensitivity_parameter,
    sensitivity_value = sensitivity_value,
    lower = lower,
    upper = upper,
    tau_hat = tau_hat,
    interval_type = interval_type,
    inference_type = inference_type,
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
    R2 = R2,
    Gamma = Gamma,
    T = T,
    R2_minus_raw = R2_minus_raw,
    R2_directional = R2_directional,
    var_full = var_full,
    var_reduced = var_reduced,
    gamma_raw = gamma_raw,
    ratio_max = ratio_max,
    ratio_q95 = ratio_q95,
    ratio_q99 = ratio_q99,
    n_ratio = n_ratio,
    l1_distance = l1_distance,
    T_control = T_control,
    T_treated = T_treated,
    T_all = T_all,
    n_T = n_T,
    full_weight_mass_T = full_weight_mass_T,
    reduced_weight_mass_T = reduced_weight_mass_T,
    max_abs_probability_shift = max_abs_probability_shift,
    q95_abs_probability_shift = q95_abs_probability_shift,
    cor_full_reduced_prob = cor_full_reduced_prob,
    stringsAsFactors = FALSE
  )
}

add_rgm_benchmark_row <- function(
    out,
    v,
    tau_hat,
    t_info,
    T,
    det_interval,
    curve,
    method_label
) {

  boot_row <- rgm_select_nearest_T_row(
    curve = curve,
    T = T
  )

  if (!is.null(boot_row) && is.finite(boot_row$lower) && is.finite(boot_row$upper)) {
    lower <- boot_row$lower
    upper <- boot_row$upper
    interval_type <- paste0("bootstrap_percentile_fixed_T_", det_interval$bound_type[1])
    inference_type <- "bootstrap"
    b_lower <- boot_row$lower
    b_upper <- boot_row$upper
    b_B <- boot_row$B_requested
    b_n <- boot_row$n_valid
    b_req <- boot_row$required_valid
    b_pool <- boot_row$bootstrap_pool
    status <- paste(t_info$benchmark_status, "bootstrap_fixed_original_T", sep = " | ")
  } else {
    lower <- det_interval$lower
    upper <- det_interval$upper
    interval_type <- paste0("deterministic_fallback_", det_interval$bound_type[1])
    inference_type <- "deterministic_fallback"
    b_lower <- NA_real_
    b_upper <- NA_real_
    b_B <- NA_real_
    b_n <- NA_real_
    b_req <- NA_real_
    b_pool <- NA_character_
    status <- paste(t_info$benchmark_status, "benchmark_bootstrap_unavailable_used_deterministic_bounds", sep = " | ")
  }

  rbind(
    out,
    make_extension_benchmark_row(
      variable = v,
      removed_terms = t_info$removed_terms,
      method = method_label,
      sensitivity_parameter = "T",
      sensitivity_value = T,
      lower = lower,
      upper = upper,
      tau_hat = tau_hat,
      interval_type = interval_type,
      inference_type = inference_type,
      deterministic_lower = det_interval$lower,
      deterministic_upper = det_interval$upper,
      bootstrap_lower = b_lower,
      bootstrap_upper = b_upper,
      bootstrap_B = b_B,
      bootstrap_n_valid = b_n,
      bootstrap_required_valid = b_req,
      bootstrap_pool = b_pool,
      benchmark_status = status,
      benchmark_group_used = t_info$benchmark_group_used,
      T = T,
      l1_distance = t_info$l1_distance,
      T_control = t_info$T_control,
      T_treated = t_info$T_treated,
      T_all = t_info$T_all,
      n_T = t_info$n_T,
      full_weight_mass_T = t_info$full_weight_mass,
      reduced_weight_mass_T = t_info$reduced_weight_mass,
      max_abs_probability_shift = t_info$max_abs_probability_shift,
      q95_abs_probability_shift = t_info$q95_abs_probability_shift,
      cor_full_reduced_prob = t_info$cor_full_reduced_prob
    )
  )
}


rgm_get_scalar <- function(row, name, default = NA) {
  if (!name %in% names(row)) {
    return(default)
  }
  value <- row[[name]][1]
  if (length(value) == 0 || is.null(value)) {
    return(default)
  }
  value
}

coerce_existing_vbm_msm_benchmark_rows <- function(base_benchmark_df) {
  if (is.null(base_benchmark_df) || nrow(base_benchmark_df) == 0) {
    return(data.frame())
  }

  keep <- base_benchmark_df$method %in% c(
    "VBM benchmark (bootstrap)",
    "MSM benchmark (bootstrap)"
  )

  base_core <- base_benchmark_df[keep, , drop = FALSE]

  if (nrow(base_core) == 0) {
    return(data.frame())
  }

  out <- data.frame()

  for (i in seq_len(nrow(base_core))) {
    row <- base_core[i, , drop = FALSE]
    raw_method <- as.character(rgm_get_scalar(row, "method", ""))

    if (identical(raw_method, "VBM benchmark (bootstrap)")) {
      method <- "VBM"
      sensitivity_parameter <- "R2"
      R2 <- suppressWarnings(as.numeric(rgm_get_scalar(row, "sensitivity_value", NA_real_)))
      Gamma <- NA_real_
    } else if (identical(raw_method, "MSM benchmark (bootstrap)")) {
      method <- "MSM"
      sensitivity_parameter <- "Gamma"
      R2 <- NA_real_
      Gamma <- suppressWarnings(as.numeric(rgm_get_scalar(row, "sensitivity_value", NA_real_)))
    } else {
      next
    }

    out <- rbind(
      out,
      make_extension_benchmark_row(
        variable = as.character(rgm_get_scalar(row, "variable", NA_character_)),
        removed_terms = as.character(rgm_get_scalar(row, "removed_terms", NA_character_)),
        method = method,
        sensitivity_parameter = sensitivity_parameter,
        sensitivity_value = suppressWarnings(as.numeric(rgm_get_scalar(row, "sensitivity_value", NA_real_))),
        lower = suppressWarnings(as.numeric(rgm_get_scalar(row, "lower", NA_real_))),
        upper = suppressWarnings(as.numeric(rgm_get_scalar(row, "upper", NA_real_))),
        tau_hat = suppressWarnings(as.numeric(rgm_get_scalar(row, "tau_hat", NA_real_))),
        interval_type = as.character(rgm_get_scalar(row, "interval_type", "bootstrap_percentile_fixed_sensitivity")),
        inference_type = as.character(rgm_get_scalar(row, "inference_type", "bootstrap")),
        deterministic_lower = NA_real_,
        deterministic_upper = NA_real_,
        bootstrap_lower = suppressWarnings(as.numeric(rgm_get_scalar(row, "lower", NA_real_))),
        bootstrap_upper = suppressWarnings(as.numeric(rgm_get_scalar(row, "upper", NA_real_))),
        bootstrap_B = suppressWarnings(as.numeric(rgm_get_scalar(row, "bootstrap_B", NA_real_))),
        bootstrap_n_valid = suppressWarnings(as.numeric(rgm_get_scalar(row, "bootstrap_n_valid", NA_real_))),
        bootstrap_required_valid = suppressWarnings(as.numeric(rgm_get_scalar(row, "bootstrap_required_valid", NA_real_))),
        bootstrap_pool = as.character(rgm_get_scalar(row, "bootstrap_pool", NA_character_)),
        benchmark_status = as.character(rgm_get_scalar(row, "benchmark_status", "reused_from_main_benchmark")),
        benchmark_group_used = as.character(rgm_get_scalar(row, "benchmark_group_used", NA_character_)),
        R2 = R2,
        Gamma = Gamma,
        R2_minus_raw = suppressWarnings(as.numeric(rgm_get_scalar(row, "R2_minus_raw", NA_real_))),
        R2_directional = suppressWarnings(as.numeric(rgm_get_scalar(row, "R2_directional", NA_real_))),
        var_full = suppressWarnings(as.numeric(rgm_get_scalar(row, "var_full", NA_real_))),
        var_reduced = suppressWarnings(as.numeric(rgm_get_scalar(row, "var_reduced", NA_real_))),
        gamma_raw = suppressWarnings(as.numeric(rgm_get_scalar(row, "gamma_raw", NA_real_))),
        ratio_max = suppressWarnings(as.numeric(rgm_get_scalar(row, "ratio_max", NA_real_))),
        ratio_q95 = suppressWarnings(as.numeric(rgm_get_scalar(row, "ratio_q95", NA_real_))),
        ratio_q99 = suppressWarnings(as.numeric(rgm_get_scalar(row, "ratio_q99", NA_real_))),
        n_ratio = suppressWarnings(as.numeric(rgm_get_scalar(row, "n_ratio", NA_real_)))
      )
    )
  }

  out
}

generate_rgm_benchmark_comparison_data <- function(
    analysis,
    tau_hat,
    config,
    full_formula = ps_formula,
    base_benchmark_df = NULL,
    base_vbm_benchmark_bootstrap_curve = NULL,
    base_msm_benchmark_bootstrap_curve = NULL
) {

  out <- coerce_existing_vbm_msm_benchmark_rows(base_benchmark_df)
  benchmark_groups <- resolve_benchmark_groups(config)

  benchmark_items <- list()
  T_values <- c()

  for (v in names(benchmark_groups)) {

    remove_terms <- benchmark_groups[[v]]

    t_info <- rgm_compute_benchmark_T(
      full_formula = full_formula,
      benchmark_name = v,
      remove_terms = remove_terms,
      data = analysis,
      config = config
    )

    T <- t_info$T_used
    T_values <- c(T_values, T)

    rgm_sharp_det <- if (is.finite(T)) {
      rgm_sharp_att_bounds(
        analysis = analysis,
        T = T,
        config = config
      )
    } else {
      data.frame(lower = NA_real_, upper = NA_real_, bound_type = "finite_sample_sharp")
    }

    rgm_conservative_det <- if (is.finite(T)) {
      rgm_conservative_att_bounds(
        analysis = analysis,
        T = T,
        config = config
      )
    } else {
      data.frame(lower = NA_real_, upper = NA_real_, bound_type = "conservative_range")
    }

    benchmark_items[[v]] <- list(
      t_info = t_info,
      T = T,
      rgm_sharp_det = rgm_sharp_det,
      rgm_conservative_det = rgm_conservative_det
    )
  }

  rgm_sharp_curve <- generate_rgm_sharp_benchmark_bootstrap_curve(
    analysis = analysis,
    T_values = T_values,
    config = config
  )

  rgm_conservative_curve <- generate_rgm_conservative_benchmark_bootstrap_curve(
    analysis = analysis,
    T_values = T_values,
    config = config
  )

  for (v in names(benchmark_items)) {

    item <- benchmark_items[[v]]

    out <- add_rgm_benchmark_row(
      out = out,
      v = v,
      tau_hat = tau_hat,
      t_info = item$t_info,
      T = item$T,
      det_interval = item$rgm_conservative_det,
      curve = rgm_conservative_curve,
      method_label = "RGM-conservative"
    )

    out <- add_rgm_benchmark_row(
      out = out,
      v = v,
      tau_hat = tau_hat,
      t_info = item$t_info,
      T = item$T,
      det_interval = item$rgm_sharp_det,
      curve = rgm_sharp_curve,
      method_label = "RGM-sharp"
    )
  }

  attr(out, "vbm_benchmark_bootstrap_curve") <- base_vbm_benchmark_bootstrap_curve
  attr(out, "msm_benchmark_bootstrap_curve") <- base_msm_benchmark_bootstrap_curve
  attr(out, "rgm_sharp_benchmark_bootstrap_curve") <- rgm_sharp_curve
  attr(out, "rgm_conservative_benchmark_bootstrap_curve") <- rgm_conservative_curve

  out
}

run_rgm_benchmark_comparison <- function(
    analysis,
    tau_hat,
    full_formula = ps_formula,
    config,
    base_benchmark_df = NULL,
    base_vbm_benchmark_bootstrap_curve = NULL,
    base_msm_benchmark_bootstrap_curve = NULL
) {

  benchmark_df <- generate_rgm_benchmark_comparison_data(
    analysis = analysis,
    tau_hat = tau_hat,
    config = config,
    full_formula = full_formula,
    base_benchmark_df = base_benchmark_df,
    base_vbm_benchmark_bootstrap_curve = base_vbm_benchmark_bootstrap_curve,
    base_msm_benchmark_bootstrap_curve = base_msm_benchmark_bootstrap_curve
  )

  list(
    benchmark_df = benchmark_df,
    vbm_benchmark_bootstrap_curve = attr(benchmark_df, "vbm_benchmark_bootstrap_curve"),
    msm_benchmark_bootstrap_curve = attr(benchmark_df, "msm_benchmark_bootstrap_curve"),
    rgm_sharp_benchmark_bootstrap_curve = attr(benchmark_df, "rgm_sharp_benchmark_bootstrap_curve"),
    rgm_conservative_benchmark_bootstrap_curve = attr(benchmark_df, "rgm_conservative_benchmark_bootstrap_curve")
  )
}

plot_rgm_benchmark_comparison <- function(
    df,
    tau_hat = NA_real_,
    config = NULL
) {

  plot_df <- df

  plot_df$method_label <- factor(
    plot_df$method,
    levels = c("MSM", "VBM", "RGM-conservative", "RGM-sharp")
  )

  if (!is.null(config$benchmark_plot_labels)) {
    label_map <- config$benchmark_plot_labels
    plot_df$variable_label <- ifelse(
      plot_df$variable %in% names(label_map),
      unname(label_map[plot_df$variable]),
      plot_df$variable
    )

    ordered_vars <- names(label_map)[names(label_map) %in% unique(plot_df$variable)]
    ordered_labels <- unname(label_map[ordered_vars])

    plot_df$variable_label <- factor(
      plot_df$variable_label,
      levels = ordered_labels
    )
  } else {
    plot_df$variable_label <- factor(
      plot_df$variable,
      levels = unique(plot_df$variable)
    )
  }

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = variable_label,
      ymin = lower,
      ymax = upper,
      color = method_label
    )
  ) +
    ggplot2::geom_linerange(
      position = ggplot2::position_dodge(width = 0.78),
      linewidth = 1.0
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)
    ) +
    ggplot2::labs(
      title = "Covariate benchmark: MSM, VBM, RGM, and RGM-sharp",
      subtitle = "RGM uses conservative range bounds; RGM-sharp uses finite-sample TV bounds.",
      x = "Benchmark covariate",
      y = "ATT interval",
      color = "Model"
    )

  if (is.finite(tau_hat)) {
    p <- p +
      ggplot2::geom_hline(
        yintercept = tau_hat,
        linetype = "dotted"
      )
  }

  p
}
