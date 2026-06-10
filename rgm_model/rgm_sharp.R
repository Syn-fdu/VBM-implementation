# ============================================================
# RGM finite-sample sharp TV/L1 calculations for the integrated main.R workflow
# ============================================================
#
# For controls, let p_i = w_i / sum_j w_j.  The sharp finite-sample RGM bound
# optimizes over all probability vectors q on the observed controls satisfying
#
#   q_i >= 0,  sum_i q_i = 1,  1/2 * sum_i |q_i - p_i| <= T.
#
# The maximal control mean is obtained by moving at most T units of probability
# mass from the lowest outcomes to the highest outcomes.  The minimal control
# mean moves mass in the opposite direction.
# ============================================================

rgm_normalize_prob <- function(w) {

  out <- rep(NA_real_, length(w))
  ok <- is.finite(w) & w > 0

  if (!any(ok)) {
    return(out)
  }

  total <- sum(w[ok])

  if (!is.finite(total) || total <= 0) {
    return(out)
  }

  out[ok] <- w[ok] / total
  out
}

rgm_tv_max_gain <- function(y, p, T) {

  T <- rgm_clamp_T(T)

  ok <- is.finite(y) & is.finite(p) & p >= 0
  y <- y[ok]
  p <- p[ok]

  if (length(y) == 0) {
    return(
      list(
        gain = NA_real_,
        moved_mass = NA_real_,
        status = "empty_support",
        n = 0L
      )
    )
  }

  p_sum <- sum(p)
  if (!is.finite(p_sum) || p_sum <= 0) {
    return(
      list(
        gain = NA_real_,
        moved_mass = NA_real_,
        status = "invalid_probability_mass",
        n = length(y)
      )
    )
  }

  p <- p / p_sum

  if (T <= 0 || length(y) == 1) {
    return(
      list(
        gain = 0,
        moved_mass = 0,
        status = "T_equals_0_or_single_support_point",
        n = length(y)
      )
    )
  }

  source_order <- order(y, decreasing = FALSE)
  sink_order <- order(y, decreasing = TRUE)

  source_cap <- p[source_order]
  sink_cap <- 1 - p[sink_order]

  source_y <- y[source_order]
  sink_y <- y[sink_order]

  i <- 1L
  j <- 1L
  moved <- 0
  gain <- 0
  tol <- 1e-12

  while (i <= length(source_order) && j <= length(sink_order) && moved < T - tol) {

    if (!is.finite(source_cap[i]) || source_cap[i] <= tol) {
      i <- i + 1L
      next
    }

    if (!is.finite(sink_cap[j]) || sink_cap[j] <= tol) {
      j <- j + 1L
      next
    }

    if (sink_y[j] <= source_y[i] + tol) {
      break
    }

    amount <- min(source_cap[i], sink_cap[j], T - moved)

    if (!is.finite(amount) || amount <= tol) {
      break
    }

    gain <- gain + amount * (sink_y[j] - source_y[i])
    moved <- moved + amount

    source_cap[i] <- source_cap[i] - amount
    sink_cap[j] <- sink_cap[j] - amount
  }

  list(
    gain = gain,
    moved_mass = moved,
    status = "sharp_greedy_transport_solution",
    n = length(y)
  )
}

rgm_sharp_control_mean_bounds <- function(y, w, T) {

  ok <- is.finite(y) & is.finite(w) & w > 0
  y <- y[ok]
  w <- w[ok]

  if (length(y) == 0) {
    return(
      list(
        control_mean_min = NA_real_,
        control_mean_max = NA_real_,
        observed_control_mean = NA_real_,
        max_gain = NA_real_,
        min_gain = NA_real_,
        moved_mass_max = NA_real_,
        moved_mass_min = NA_real_,
        n_control = 0L,
        status = "empty_control_set"
      )
    )
  }

  p <- rgm_normalize_prob(w)
  ok_p <- is.finite(p) & p >= 0
  y <- y[ok_p]
  p <- p[ok_p]

  if (length(y) == 0 || sum(p) <= 0) {
    return(
      list(
        control_mean_min = NA_real_,
        control_mean_max = NA_real_,
        observed_control_mean = NA_real_,
        max_gain = NA_real_,
        min_gain = NA_real_,
        moved_mass_max = NA_real_,
        moved_mass_min = NA_real_,
        n_control = 0L,
        status = "invalid_probability_weights"
      )
    )
  }

  p <- p / sum(p)
  observed_mean <- sum(p * y)

  max_info <- rgm_tv_max_gain(
    y = y,
    p = p,
    T = T
  )

  min_info <- rgm_tv_max_gain(
    y = -y,
    p = p,
    T = T
  )

  control_mean_max <- observed_mean + max_info$gain
  control_mean_min <- observed_mean - min_info$gain

  list(
    control_mean_min = control_mean_min,
    control_mean_max = control_mean_max,
    observed_control_mean = observed_mean,
    max_gain = max_info$gain,
    min_gain = min_info$gain,
    moved_mass_max = max_info$moved_mass,
    moved_mass_min = min_info$moved_mass,
    n_control = length(y),
    status = paste(max_info$status, min_info$status, sep = " | ")
  )
}

rgm_sharp_att_bounds <- function(analysis, T, config = NULL) {

  if (!all(c("Y", "Z", "w") %in% names(analysis))) {
    stop("analysis must contain Y, Z, and w columns.")
  }

  T <- rgm_clamp_T(T)

  Y_treated <- analysis$Y[analysis$Z == 1]
  Y0 <- analysis$Y[analysis$Z == 0]
  w0 <- analysis$w[analysis$Z == 0]

  Y_treated <- Y_treated[is.finite(Y_treated)]

  if (length(Y_treated) == 0) {
    return(
      data.frame(
        method = "RGM-sharp",
        bound_type = "finite_sample_sharp",
        T = T,
        lower = NA_real_,
        upper = NA_real_,
        tau_hat = NA_real_,
        treated_mean = NA_real_,
        observed_control_mean = NA_real_,
        control_mean_min = NA_real_,
        control_mean_max = NA_real_,
        max_gain = NA_real_,
        min_gain = NA_real_,
        moved_mass_max = NA_real_,
        moved_mass_min = NA_real_,
        n_treated = 0L,
        n_control = sum(analysis$Z == 0, na.rm = TRUE),
        status = "empty_treated_set",
        stringsAsFactors = FALSE
      )
    )
  }

  treated_mean <- mean(Y_treated)

  cm <- rgm_sharp_control_mean_bounds(
    y = Y0,
    w = w0,
    T = T
  )

  tau_hat <- treated_mean - cm$observed_control_mean

  data.frame(
    method = "RGM-sharp",
    bound_type = "finite_sample_sharp",
    T = T,
    lower = treated_mean - cm$control_mean_max,
    upper = treated_mean - cm$control_mean_min,
    tau_hat = tau_hat,
    treated_mean = treated_mean,
    observed_control_mean = cm$observed_control_mean,
    control_mean_min = cm$control_mean_min,
    control_mean_max = cm$control_mean_max,
    max_gain = cm$max_gain,
    min_gain = cm$min_gain,
    moved_mass_max = cm$moved_mass_max,
    moved_mass_min = cm$moved_mass_min,
    n_treated = length(Y_treated),
    n_control = cm$n_control,
    status = cm$status,
    stringsAsFactors = FALSE
  )
}

run_rgm_sharp_att_T_grid <- function(analysis, T_grid, config = NULL) {

  out <- data.frame()

  for (t in T_grid) {
    bb <- rgm_sharp_att_bounds(
      analysis = analysis,
      T = t,
      config = config
    )
    out <- rbind(out, bb)
  }

  out
}

# ------------------------------------------------------------
# Sharp bootstrap ATT-T curve and T* search
# ------------------------------------------------------------

generate_rgm_sharp_bootstrap_curve <- function(
    data,
    T_grid,
    config,
    B = NULL,
    stage = "sharp"
) {

  if (is.null(B)) {
    B <- rgm_ext_value(
      config,
      "rgm_sharp_bootstrap_B",
      rgm_ext_value(config, "fine_B", 300)
    )
  }

  if (is.null(B) || !is.finite(B) || B <= 0) {
    return(data.frame())
  }

  B <- as.integer(B)
  T_grid <- rgm_clamp_T(T_grid)

  lower_mat <- matrix(NA_real_, nrow = B, ncol = length(T_grid))
  upper_mat <- matrix(NA_real_, nrow = B, ncol = length(T_grid))

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

    for (j in seq_along(T_grid)) {
      bb <- tryCatch(
        rgm_sharp_att_bounds(
          analysis = d,
          T = T_grid[j],
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

  for (j in seq_along(T_grid)) {

    lower_vals <- lower_mat[, j]
    upper_vals <- upper_mat[, j]

    ok <- is.finite(lower_vals) & is.finite(upper_vals)
    n_valid <- sum(ok)

    if (n_valid >= required_valid) {
      lower <- as.numeric(
        stats::quantile(
          lower_vals[ok],
          probs = config$percentile[1],
          na.rm = TRUE,
          names = FALSE
        )
      )

      upper <- as.numeric(
        stats::quantile(
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
        method = "RGM-sharp",
        bound_type = "finite_sample_sharp",
        stage = stage,
        T = T_grid[j],
        lower = lower,
        upper = upper,
        n_valid = n_valid,
        B_requested = B,
        required_valid = required_valid,
        bootstrap_pool = paste0("rgm_sharp_stage_", stage),
        stringsAsFactors = FALSE
      )
    )
  }

  attr(out, "failures") <- failures
  attr(out, "B_requested") <- B
  out
}

rgm_interpolate_T_crossing <- function(results) {

  results <- results[order(results$T), , drop = FALSE]

  crossing_index <- which(
    is.finite(results$lower) &
      results$lower <= 0
  )[1]

  if (is.na(crossing_index)) {
    return(NA_real_)
  }

  if (crossing_index == 1) {
    return(results$T[crossing_index])
  }

  previous_candidates <- which(
    is.finite(results$lower) &
      seq_len(nrow(results)) < crossing_index
  )

  if (length(previous_candidates) == 0) {
    return(results$T[crossing_index])
  }

  previous <- results[tail(previous_candidates, 1), , drop = FALSE]
  current <- results[crossing_index, , drop = FALSE]

  if (!is.finite(previous$lower) || !is.finite(current$lower)) {
    return(current$T)
  }

  if (previous$lower <= 0) {
    return(current$T)
  }

  previous$T +
    (0 - previous$lower) /
    (current$lower - previous$lower) *
    (current$T - previous$T)
}

find_rgm_sharp_T_bootstrap <- function(data, config) {

  seed_offset <- rgm_ext_value(config, "rgm_sharp_search_seed_offset", 12100)
  if (!is.null(config$seed) && is.finite(config$seed)) {
    set.seed(as.integer(config$seed + seed_offset))
  }

  T_min <- rgm_ext_value(config, "rgm_T_min", rgm_ext_value(config, "T_min", 0))
  T_max <- rgm_ext_value(config, "rgm_T_max", rgm_ext_value(config, "T_max", 1))
  T_coarse_step <- rgm_ext_value(config, "rgm_T_coarse_step", rgm_ext_value(config, "T_coarse_step", 0.05))
  T_fine_step <- rgm_ext_value(config, "rgm_T_fine_step", rgm_ext_value(config, "T_fine_step", 0.005))
  T_fine_window <- rgm_ext_value(config, "rgm_T_fine_window", rgm_ext_value(config, "T_fine_window", 0.05))

  coarse_B <- rgm_ext_value(
    config,
    "rgm_sharp_coarse_B",
    rgm_ext_value(config, "coarse_B", 500)
  )

  fine_B <- rgm_ext_value(
    config,
    "rgm_sharp_fine_B",
    rgm_ext_value(config, "fine_B", 1000)
  )

  coarse_grid <- seq(
    T_min,
    T_max,
    by = T_coarse_step
  )

  cat(
    "Using sharp RGM bootstrap pools: coarse B =",
    coarse_B,
    ", fine B =",
    fine_B,
    "\n"
  )

  coarse_results <- generate_rgm_sharp_bootstrap_curve(
    data = data,
    T_grid = coarse_grid,
    config = config,
    B = as.integer(coarse_B),
    stage = "coarse"
  )

  coarse_crossing <- which(
    is.finite(coarse_results$lower) &
      coarse_results$lower <= 0
  )[1]

  diagnostics <- data.frame(
    stage = "coarse",
    B_requested = as.integer(coarse_B),
    failure_one_class = as.integer(attr(coarse_results, "failures")["one_class"]),
    failure_ps_fit = as.integer(attr(coarse_results, "failures")["ps_fit"]),
    failure_weights = as.integer(attr(coarse_results, "failures")["weights"]),
    failure_bounds = as.integer(attr(coarse_results, "failures")["bounds"]),
    stringsAsFactors = FALSE
  )

  if (is.na(coarse_crossing)) {
    return(
      list(
        results = coarse_results,
        T_star = NA_real_,
        coarse_results = coarse_results,
        fine_results = NULL,
        diagnostics = diagnostics
      )
    )
  }

  coarse_t <- coarse_results$T[coarse_crossing]

  fine_start <- max(T_min, coarse_t - T_fine_window)
  fine_end <- min(T_max, coarse_t + T_fine_window)

  fine_grid <- seq(
    fine_start,
    fine_end,
    by = T_fine_step
  )

  fine_results <- generate_rgm_sharp_bootstrap_curve(
    data = data,
    T_grid = fine_grid,
    config = config,
    B = as.integer(fine_B),
    stage = "fine"
  )

  diagnostics <- rbind(
    diagnostics,
    data.frame(
      stage = "fine",
      B_requested = as.integer(fine_B),
      failure_one_class = as.integer(attr(fine_results, "failures")["one_class"]),
      failure_ps_fit = as.integer(attr(fine_results, "failures")["ps_fit"]),
      failure_weights = as.integer(attr(fine_results, "failures")["weights"]),
      failure_bounds = as.integer(attr(fine_results, "failures")["bounds"]),
      stringsAsFactors = FALSE
    )
  )

  T_star <- rgm_interpolate_T_crossing(fine_results)

  list(
    results = rbind(coarse_results, fine_results),
    T_star = T_star,
    coarse_results = coarse_results,
    fine_results = fine_results,
    diagnostics = diagnostics
  )
}

generate_rgm_sharp_benchmark_bootstrap_curve <- function(
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
    "rgm_sharp_benchmark_seed_offset",
    13100
  )

  if (!is.null(config$seed) && is.finite(config$seed)) {
    set.seed(as.integer(config$seed + seed_offset))
  }

  curve <- generate_rgm_sharp_bootstrap_curve(
    data = analysis,
    T_grid = T_values,
    config = config,
    B = as.integer(B),
    stage = "benchmark"
  )

  if (nrow(curve) > 0) {
    curve$bootstrap_pool <- "shared_fixed_benchmark_T_pool_sharp"
    curve$benchmark_bootstrap_type <- "fixed_original_sample_T_sharp"
  }

  curve
}

plot_rgm_sharp_att_T_curve <- function(
    results,
    T_star = NA_real_,
    deterministic_curve = NULL
) {

  p <- ggplot2::ggplot(results, ggplot2::aes(x = T)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      fill = "grey60",
      alpha = 0.5
    ) +
    ggplot2::geom_line(ggplot2::aes(y = lower), linewidth = 0.55) +
    ggplot2::geom_line(ggplot2::aes(y = upper), linewidth = 0.55) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "red",
      linetype = "dashed"
    ) +
    ggplot2::facet_wrap(~stage) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "RGM sharp finite-sample bootstrap ATT-T curve",
      x = "T (total-variation/L1 distance)",
      y = "ATT interval"
    )

  if (!is.null(deterministic_curve) && nrow(deterministic_curve) > 0) {
    p <- p +
      ggplot2::geom_line(
        data = deterministic_curve,
        ggplot2::aes(x = T, y = lower),
        inherit.aes = FALSE,
        linetype = "dotted"
      ) +
      ggplot2::geom_line(
        data = deterministic_curve,
        ggplot2::aes(x = T, y = upper),
        inherit.aes = FALSE,
        linetype = "dotted"
      )
  }

  if (is.finite(T_star)) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = T_star,
        linetype = "dotted",
        color = "red"
      )
  }

  p
}
