# ============================================================
# Bootstrap variance-based sensitivity analysis
# ============================================================
#
# Bootstrap re-sampling and propensity-score re-fitting are done once per
# search stage, not once per R2 value.  The resulting {tau_b, scale_b, base_scale_b} pool
# is reused across all R2 values in that stage, which keeps the curve smooth
# without changing the percentile-bootstrap target.
# ============================================================

.config_value <- function(config, name, default) {
  if (is.null(config) || is.null(config[[name]])) {
    return(default)
  }
  config[[name]]
}

bootstrap_resample_data <- function(data, config) {

  resample_by_stage <- isTRUE(
    .config_value(config, "bootstrap_resample_by_stage", TRUE)
  )

  if (resample_by_stage && "Z" %in% names(data)) {

    idx1 <- which(data$Z == 1)
    idx0 <- which(data$Z == 0)

    if (length(idx1) > 0 && length(idx0) > 0) {
      idx <- c(
        sample(idx1, size = length(idx1), replace = TRUE),
        sample(idx0, size = length(idx0), replace = TRUE)
      )

      return(data[idx, , drop = FALSE])
    }
  }

  idx <- sample.int(
    n = nrow(data),
    size = nrow(data),
    replace = TRUE
  )

  data[idx, , drop = FALSE]
}

bootstrap_required_valid <- function(B, config) {

  min_fraction <- .config_value(
    config,
    "min_bootstrap_valid_fraction",
    0.5
  )

  min_abs <- .config_value(
    config,
    "min_bootstrap_valid",
    20
  )

  min(
    B,
    max(
      min_abs,
      ceiling(min_fraction * B)
    )
  )
}

generate_bootstrap_stats <- function(
    data,
    B = 100,
    config
) {

  tau_vals <- rep(NA_real_, B)
  scale_vals <- rep(NA_real_, B)
  base_scale_vals <- rep(NA_real_, B)

  failures <- c(
    one_class = 0,
    ps_fit = 0,
    weights = 0,
    tau_or_scale = 0
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

    scale_b <- bias_scale(
      w = w_b,
      Y = d$Y,
      Z = d$Z,
      config = config
    )

    base_scale_b <- bias_scale_base(
      w = w_b,
      Y = d$Y,
      Z = d$Z,
      config = config
    )

    if (!is.finite(tau_b) || !is.finite(scale_b) || !is.finite(base_scale_b)) {
      failures["tau_or_scale"] <- failures["tau_or_scale"] + 1
      next
    }

    tau_vals[b] <- tau_b
    scale_vals[b] <- scale_b
    base_scale_vals[b] <- base_scale_b
  }

  stats <- data.frame(
    b = seq_len(B),
    tau = tau_vals,
    scale = scale_vals,
    base_scale = base_scale_vals
  )

  stats <- stats[
    is.finite(stats$tau) &
      is.finite(stats$scale) &
      is.finite(stats$base_scale),
    ,
    drop = FALSE
  ]

  attr(stats, "B_requested") <- B
  attr(stats, "failures") <- failures

  stats
}

bootstrap_ci_from_stats <- function(
    stats,
    R2,
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

  lower_vals <- stats$tau - stats$scale * multiplier
  upper_vals <- stats$tau + stats$scale * multiplier

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

run_bootstrap_grid_from_stats <- function(
    stats,
    R2_grid,
    stage,
    config,
    pool_label = NULL
) {

  if (is.null(pool_label)) {
    pool_label <- paste0("stage_specific_", stage)
  }

  B_effective <- nrow(stats)
  results <- data.frame()

  for (r in R2_grid) {

    ci <- bootstrap_ci_from_stats(
      stats = stats,
      R2 = r,
      B = B_effective,
      config = config
    )

    results <- rbind(
      results,
      data.frame(
        stage = stage,
        R2 = r,
        lower = ci["lower"],
        upper = ci["upper"],
        n_valid = ci["n_valid"],
        B = B_effective,
        bootstrap_pool = pool_label,
        stringsAsFactors = FALSE
      )
    )

    cat(
      stage,
      "R² =", round(r, 3),
      "CI = [",
      round(ci["lower"], 3),
      ",",
      round(ci["upper"], 3),
      "]",
      "valid =", ci["n_valid"],
      "| bootstrap pool:", pool_label, "\n"
    )
  }

  results
}

run_bootstrap_grid_stage_specific <- function(
    data,
    R2_grid,
    B,
    stage,
    config
) {

  stats <- generate_bootstrap_stats(
    data = data,
    B = B,
    config = config
  )

  if (isTRUE(.config_value(config, "diagnostics", FALSE))) {
    print_bootstrap_stats_message(stats, B, stage)
  }

  pool_label <- paste0("stage_specific_", stage)

  results <- run_bootstrap_grid_from_stats(
    stats = stats,
    R2_grid = R2_grid,
    stage = stage,
    config = config,
    pool_label = pool_label
  )

  attr(results, "bootstrap_stats_summary") <- bootstrap_stats_summary(
    stats = stats,
    stage = stage,
    scope = pool_label
  )

  results
}

print_bootstrap_stats_message <- function(stats, B, stage) {

  failures <- attr(stats, "failures")

  if (is.null(failures)) {
    failures <- c()
  }

  failure_msg <- paste(
    names(failures),
    as.integer(failures),
    sep = "=",
    collapse = ", "
  )

  cat(
    stage,
    "bootstrap stats:",
    nrow(stats),
    "valid out of",
    B,
    "| failures:",
    failure_msg,
    "\n"
  )
}

bootstrap_stats_summary <- function(
    stats,
    stage,
    scope
) {

  failures <- attr(stats, "failures")
  B_requested <- attr(stats, "B_requested")

  if (is.null(failures)) {
    failures <- c(
      one_class = NA_integer_,
      ps_fit = NA_integer_,
      weights = NA_integer_,
      tau_or_scale = NA_integer_
    )
  }

  data.frame(
    stage = stage,
    scope = scope,
    B_requested = B_requested,
    n_valid = nrow(stats),
    failure_one_class = as.integer(failures["one_class"]),
    failure_ps_fit = as.integer(failures["ps_fit"]),
    failure_weights = as.integer(failures["weights"]),
    failure_tau_or_scale = as.integer(failures["tau_or_scale"]),
    tau_mean = mean(stats$tau, na.rm = TRUE),
    tau_sd = stats::sd(stats$tau, na.rm = TRUE),
    scale_mean = mean(stats$scale, na.rm = TRUE),
    scale_sd = stats::sd(stats$scale, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

interpolate_crossing <- function(results) {

  results <- results[order(results$R2), ]

  crossing_index <- which(
    is.finite(results$lower) &
      results$lower <= 0
  )[1]

  if (is.na(crossing_index)) {
    return(NA_real_)
  }

  if (crossing_index == 1) {
    return(results$R2[crossing_index])
  }

  previous_candidates <- which(
    is.finite(results$lower) &
      seq_len(nrow(results)) < crossing_index
  )

  if (length(previous_candidates) == 0) {
    return(results$R2[crossing_index])
  }

  previous <- results[tail(previous_candidates, 1), ]
  current <- results[crossing_index, ]

  if (!is.finite(previous$lower) || !is.finite(current$lower)) {
    return(current$R2)
  }

  if (previous$lower <= 0) {
    return(current$R2)
  }

  previous$R2 +
    (0 - previous$lower) /
    (current$lower - previous$lower) *
    (current$R2 - previous$R2)
}

find_R2_bootstrap <- function(
    data,
    config
) {

  set.seed(config$seed)

  coarse_grid <- seq(
    config$R2_min,
    config$R2_max,
    by = config$R2_coarse_step
  )

  cat(
    "Using stage-specific shared bootstrap pools: coarse B =",
    config$coarse_B,
    ", fine B =",
    config$fine_B,
    "\n"
  )

  coarse_results <- run_bootstrap_grid_stage_specific(
    data = data,
    R2_grid = coarse_grid,
    B = config$coarse_B,
    stage = "coarse",
    config = config
  )

  stats_summary <- attr(coarse_results, "bootstrap_stats_summary")

  coarse_crossing <- which(
    is.finite(coarse_results$lower) &
      coarse_results$lower <= 0
  )[1]

  if (is.na(coarse_crossing)) {
    return(
      list(
        results = coarse_results,
        R2_star = NA_real_,
        coarse_results = coarse_results,
        fine_results = NULL,
        bootstrap_stats_summary = stats_summary
      )
    )
  }

  coarse_r <- coarse_results$R2[coarse_crossing]

  fine_start <- max(
    config$R2_min,
    coarse_r - config$R2_fine_window
  )

  fine_end <- min(
    config$R2_max,
    coarse_r + config$R2_fine_window
  )

  fine_grid <- seq(
    fine_start,
    fine_end,
    by = config$R2_fine_step
  )

  fine_results <- run_bootstrap_grid_stage_specific(
    data = data,
    R2_grid = fine_grid,
    B = config$fine_B,
    stage = "fine",
    config = config
  )

  stats_summary <- rbind(
    stats_summary,
    attr(fine_results, "bootstrap_stats_summary")
  )

  R2_star <- interpolate_crossing(fine_results)

  list(
    results = rbind(coarse_results, fine_results),
    R2_star = R2_star,
    coarse_results = coarse_results,
    fine_results = fine_results,
    bootstrap_stats_summary = stats_summary
  )
}
