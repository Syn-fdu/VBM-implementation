# ============================================================
# Marginal Sensitivity Model (MSM) for ATT weighting estimators
# ============================================================
#
# This file implements the standard marginal sensitivity model for the
# control weights in an inverse-propensity-weighted ATT estimator:
#
#   Gamma^{-1} <= lambda_i = w_i^* / w_i <= Gamma.
#
# For a fixed Gamma, the worst-case imputed control mean is the sharp
# finite-sample solution to
#
#   max/min_a  sum_i a_i w_i Y_i / sum_i a_i w_i
#   subject to a_i in [1/Gamma, Gamma].
#
# The optimum is obtained by assigning the largest multiplier to one tail
# of the control outcomes and the smallest multiplier to the other tail.
# ============================================================

normalize_msm_weights <- function(w) {

  w <- w[is.finite(w) & w > 0]

  if (length(w) == 0) {
    return(w)
  }

  # Scale normalization only improves numerical conditioning; the ratio
  # objective itself is invariant to multiplying all weights by a constant.
  w / mean(w)
}

msm_sharp_weighted_mean <- function(
    y,
    w,
    Gamma,
    direction = c("max", "min")
) {

  direction <- match.arg(direction)

  if (!is.finite(Gamma) || Gamma < 1) {
    stop("Gamma must be finite and >= 1.")
  }

  ok <- is.finite(y) & is.finite(w) & w > 0
  y <- y[ok]
  w <- w[ok]

  if (length(y) == 0) {
    return(
      list(
        mean = NA_real_,
        threshold_index = NA_integer_,
        n_control = 0,
        status = "empty_control_set"
      )
    )
  }

  w <- normalize_msm_weights(w)

  if (Gamma == 1) {
    return(
      list(
        mean = stats::weighted.mean(y, w),
        threshold_index = NA_integer_,
        n_control = length(y),
        status = "Gamma_equals_1_observed_weighted_mean"
      )
    )
  }

  if (direction == "max") {
    ord <- order(y, decreasing = TRUE)
  } else {
    ord <- order(y, decreasing = FALSE)
  }

  y_ord <- y[ord]
  w_ord <- w[ord]

  lo <- 1 / Gamma
  hi <- Gamma

  cum_wy <- c(0, cumsum(w_ord * y_ord))
  cum_w <- c(0, cumsum(w_ord))

  total_wy <- sum(w_ord * y_ord)
  total_w <- sum(w_ord)

  numerator <- hi * cum_wy + lo * (total_wy - cum_wy)
  denominator <- hi * cum_w + lo * (total_w - cum_w)

  values <- numerator / denominator
  values[!is.finite(values)] <- NA_real_

  if (all(is.na(values))) {
    return(
      list(
        mean = NA_real_,
        threshold_index = NA_integer_,
        n_control = length(y),
        status = "all_threshold_values_invalid"
      )
    )
  }

  if (direction == "max") {
    best <- which.max(values)
  } else {
    best <- which.min(values)
  }

  list(
    mean = as.numeric(values[best]),
    threshold_index = as.integer(best - 1),
    n_control = length(y),
    status = "sharp_fractional_linear_solution"
  )
}

msm_att_bounds <- function(
    analysis,
    Gamma,
    config = NULL
) {

  if (!("Y" %in% names(analysis)) || !("Z" %in% names(analysis)) || !("w" %in% names(analysis))) {
    stop("analysis must contain Y, Z, and w columns.")
  }

  Y_treated <- analysis$Y[analysis$Z == 1]
  Y0 <- analysis$Y[analysis$Z == 0]
  w0 <- analysis$w[analysis$Z == 0]

  ok_t <- is.finite(Y_treated)
  Y_treated <- Y_treated[ok_t]

  ok0 <- is.finite(Y0) & is.finite(w0) & w0 > 0
  Y0 <- Y0[ok0]
  w0 <- w0[ok0]

  if (length(Y_treated) == 0 || length(Y0) == 0) {
    return(
      data.frame(
        Gamma = Gamma,
        lower = NA_real_,
        upper = NA_real_,
        control_mean_min = NA_real_,
        control_mean_max = NA_real_,
        treated_mean = NA_real_,
        n_treated = length(Y_treated),
        n_control = length(Y0),
        status = "empty_treated_or_control",
        stringsAsFactors = FALSE
      )
    )
  }

  if (!is.null(config$msm_min_control_n) &&
      is.finite(config$msm_min_control_n) &&
      length(Y0) < config$msm_min_control_n) {
    return(
      data.frame(
        Gamma = Gamma,
        lower = NA_real_,
        upper = NA_real_,
        control_mean_min = NA_real_,
        control_mean_max = NA_real_,
        treated_mean = mean(Y_treated),
        n_treated = length(Y_treated),
        n_control = length(Y0),
        status = "too_few_controls",
        stringsAsFactors = FALSE
      )
    )
  }

  treated_mean <- mean(Y_treated)

  control_min <- msm_sharp_weighted_mean(
    y = Y0,
    w = w0,
    Gamma = Gamma,
    direction = "min"
  )

  control_max <- msm_sharp_weighted_mean(
    y = Y0,
    w = w0,
    Gamma = Gamma,
    direction = "max"
  )

  data.frame(
    Gamma = Gamma,
    lower = treated_mean - control_max$mean,
    upper = treated_mean - control_min$mean,
    control_mean_min = control_min$mean,
    control_mean_max = control_max$mean,
    treated_mean = treated_mean,
    n_treated = length(Y_treated),
    n_control = length(Y0),
    min_threshold_index = control_min$threshold_index,
    max_threshold_index = control_max$threshold_index,
    status = ifelse(
      identical(control_min$status, "sharp_fractional_linear_solution") ||
        identical(control_min$status, "Gamma_equals_1_observed_weighted_mean"),
      control_max$status,
      paste(control_min$status, control_max$status, sep = " | ")
    ),
    stringsAsFactors = FALSE
  )
}

generate_msm_bootstrap_curve <- function(
    data,
    config,
    B = NULL,
    Gamma_grid = NULL
) {

  if (is.null(B)) {
    B <- config$msm_bootstrap_B
  }

  if (is.null(B) || !is.finite(B) || B <= 0) {
    return(data.frame())
  }

  B <- as.integer(B)

  if (is.null(Gamma_grid)) {
    Gamma_grid <- seq(
      config$Gamma_min,
      config$Gamma_max,
      by = config$Gamma_step
    )
  }

  lower_mat <- matrix(NA_real_, nrow = B, ncol = length(Gamma_grid))
  upper_mat <- matrix(NA_real_, nrow = B, ncol = length(Gamma_grid))

  failures <- c(
    one_class = 0,
    ps_fit = 0,
    weights = 0,
    bounds = 0
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

    d$ps <- ps_b
    d$w <- w_b

    for (j in seq_along(Gamma_grid)) {
      bb <- tryCatch(
        msm_att_bounds(
          analysis = d,
          Gamma = Gamma_grid[j],
          config = config
        ),
        error = function(e) NULL
      )

      if (is.null(bb) || !is.finite(bb$lower) || !is.finite(bb$upper)) {
        failures["bounds"] <- failures["bounds"] + 1
        next
      }

      lower_mat[b, j] <- bb$lower
      upper_mat[b, j] <- bb$upper
    }
  }

  out <- data.frame()

  required_valid <- bootstrap_required_valid(
    B = B,
    config = config
  )

  for (j in seq_along(Gamma_grid)) {

    lower_vals <- lower_mat[, j]
    upper_vals <- upper_mat[, j]

    ok <- is.finite(lower_vals) & is.finite(upper_vals)
    n_valid <- sum(ok)

    if (n_valid >= required_valid) {
      lower <- as.numeric(
        quantile(
          lower_vals[ok],
          probs = config$percentile[1],
          na.rm = TRUE,
          names = FALSE
        )
      )

      upper <- as.numeric(
        quantile(
          upper_vals[ok],
          probs = config$percentile[2],
          na.rm = TRUE,
          names = FALSE
        )
      )
    } else {
      lower <- NA_real_
      upper <- NA_real_
    }

    out <- rbind(
      out,
      data.frame(
        Gamma = Gamma_grid[j],
        lower = lower,
        upper = upper,
        n_valid = n_valid,
        B_requested = B,
        required_valid = required_valid,
        stringsAsFactors = FALSE
      )
    )
  }

  attr(out, "failures") <- failures
  out
}

make_benchmark_config <- function(config) {

  cfg <- config

  if (!is.null(config$benchmark_ps_clip_lower)) {
    cfg$ps_clip_lower <- config$benchmark_ps_clip_lower
  }

  if (!is.null(config$benchmark_ps_clip_upper)) {
    cfg$ps_clip_upper <- config$benchmark_ps_clip_upper
  }

  if (!is.null(config$benchmark_weight_truncation)) {
    cfg$weight_truncation <- config$benchmark_weight_truncation
  }

  cfg
}

msm_gamma_group_index <- function(data, config) {

  group <- config$msm_gamma_group
  if (is.null(group)) {
    group <- "control_weight_function"
  }

  if (group %in% c("treated", "treated_weight_function", "A1")) {
    return(list(idx = data$Z == 1, group = "treated_weight_function"))
  }

  if (group %in% c("control", "control_weight_function", "Z0")) {
    return(list(idx = data$Z == 0, group = "control_weight_function"))
  }

  if (group %in% c("all", "all_weight_function")) {
    return(list(idx = rep(TRUE, nrow(data)), group = "all_weight_function"))
  }

  stop(paste0("Unknown MSM Gamma group: ", group))
}

compute_benchmark_Gamma <- function(
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

  gidx <- msm_gamma_group_index(data, config)
  idx <- gidx$idx

  ratio <- w_full[idx] / w_reduced[idx]
  ratio <- ratio[is.finite(ratio) & ratio > 0]

  if (length(ratio) == 0) {
    return(
      list(
        gamma_raw = NA_real_,
        gamma_used = NA_real_,
        gamma_group_used = gidx$group,
        ratio_max = NA_real_,
        ratio_q95 = NA_real_,
        ratio_q99 = NA_real_,
        ratio_qbal = NA_real_,
        status = "empty_ratio",
        n_ratio = 0
      )
    )
  }

  log_abs_ratio <- abs(log(ratio))
  log_abs_ratio <- log_abs_ratio[is.finite(log_abs_ratio)]

  if (length(log_abs_ratio) == 0) {
    return(
      list(
        gamma_raw = NA_real_,
        gamma_used = NA_real_,
        gamma_group_used = gidx$group,
        ratio_max = NA_real_,
        ratio_q95 = NA_real_,
        ratio_q99 = NA_real_,
        ratio_qbal = NA_real_,
        status = "empty_log_ratio",
        n_ratio = 0
      )
    )
  }

  gamma_raw <- exp(max(log_abs_ratio, na.rm = TRUE))
  qbal_summary <- msm_qbal_ratio_summary(
    log_abs_ratio = log_abs_ratio,
    config = config
  )

  list(
    gamma_raw = gamma_raw,
    gamma_used = gamma_raw,
    gamma_group_used = gidx$group,
    ratio_max = gamma_raw,
    ratio_q95 = qbal_summary$ratio_q95,
    ratio_q99 = qbal_summary$ratio_q99,
    ratio_qbal = qbal_summary$ratio_qbal,
    status = "raw_max_symmetric_ratio",
    n_ratio = length(log_abs_ratio)
  )
}
