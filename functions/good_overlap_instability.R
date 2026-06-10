# ============================================================
# Current project extension utilities: instability under good overlap
# ============================================================
# Additive only: these functions do not alter the main VBM/MSM/RGM workflow.

good_overlap_sigmoid <- function(x) {
  1 / (1 + exp(-pmax(pmin(x, 35), -35)))
}

good_overlap_safe_var <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  v <- stats::var(x)
  if (!is.finite(v)) return(NA_real_)
  v
}

make_good_overlap_instability_data <- function(
    n = 1800,
    observed_strength = 0,
    hidden_strength = 0.12,
    outcome_hidden_effect = 0.75,
    treatment_effect = 0,
    seed = 90210
) {
  set.seed(seed)
  X <- stats::rnorm(n)
  U <- stats::rnorm(n)

  ps_observed_true <- good_overlap_sigmoid(observed_strength * X)
  ps_ideal_true <- good_overlap_sigmoid(observed_strength * X + hidden_strength * U)
  Z <- stats::rbinom(n, size = 1, prob = ps_ideal_true)

  Y0 <- 1 + 0.35 * X + outcome_hidden_effect * U + stats::rnorm(n)
  Y1 <- Y0 + treatment_effect
  Y <- ifelse(Z == 1, Y1, Y0)

  p1 <- min(max(mean(Z == 1), 1e-6), 1 - 1e-6)
  w_observed_true <- (ps_observed_true / (1 - ps_observed_true)) * ((1 - p1) / p1)
  w_ideal_true <- (ps_ideal_true / (1 - ps_ideal_true)) * ((1 - p1) / p1)

  out <- data.frame(
    Y = Y,
    Z = Z,
    X = X,
    U = U,
    ps_observed_true = ps_observed_true,
    ps_ideal_true = ps_ideal_true,
    w_observed_true = w_observed_true,
    w_ideal_true = w_ideal_true
  )
  attr(out, "good_overlap_params") <- list(
    n = n,
    observed_strength = observed_strength,
    hidden_strength = hidden_strength,
    outcome_hidden_effect = outcome_hidden_effect,
    treatment_effect = treatment_effect,
    seed = seed
  )
  out
}

good_overlap_metric_row <- function(
    data,
    observed_strength,
    hidden_strength,
    epsilon_fraction = 1e-4,
    variance_floor = 1e-10,
    ratio_warning_threshold = 100,
    ps_iqr_warning_threshold = 0.02,
    weight_cv_warning_threshold = 0.02
) {
  idx0 <- data$Z == 0
  w_obs0 <- data$w_observed_true[idx0]
  w_ideal0 <- data$w_ideal_true[idx0]
  y0 <- data$Y[idx0]

  var_observed <- good_overlap_safe_var(w_obs0)
  var_ideal <- good_overlap_safe_var(w_ideal0)
  mean_w <- mean(w_obs0, na.rm = TRUE)
  eps <- max(variance_floor, epsilon_fraction * mean_w^2)

  original_ratio <- ifelse(
    is.finite(var_observed) && var_observed > 0,
    var_ideal / var_observed,
    Inf
  )
  stabilized_ratio <- ifelse(
    is.finite(var_ideal),
    var_ideal / (var_observed + eps),
    NA_real_
  )

  original_R2 <- ifelse(
    is.finite(var_observed) && is.finite(var_ideal) && var_ideal > 0,
    max(0, min(1 - 1e-12, 1 - var_observed / var_ideal)),
    NA_real_
  )
  stabilized_R2 <- ifelse(
    is.finite(stabilized_ratio) && stabilized_ratio > 1,
    max(0, min(1 - 1e-12, 1 - 1 / stabilized_ratio)),
    0
  )

  cor_limit <- 1
  if (is.finite(var_observed) && var_observed > variance_floor) {
    rho <- suppressWarnings(stats::cor(w_obs0, y0))
    if (is.finite(rho)) cor_limit <- sqrt(max(0, 1 - max(-1, min(1, rho))^2))
  }

  y_var <- good_overlap_safe_var(y0)

  original_bound <- NA_real_
  if (is.finite(original_R2) && is.finite(y_var) &&
      is.finite(var_observed) && var_observed > variance_floor &&
      original_R2 > 0) {
    original_bound <- sqrt(y_var * var_observed) *
      cor_limit * sqrt(original_R2 / (1 - original_R2))
  }

  stabilized_bound <- NA_real_
  if (is.finite(stabilized_R2) && is.finite(y_var) && stabilized_R2 > 0) {
    stabilized_bound <- sqrt(y_var * (var_observed + eps)) *
      cor_limit * sqrt(stabilized_R2 / (1 - stabilized_R2))
  }

  tau_observed <- if (exists("estimate_att", mode = "function")) {
    estimate_att(data$Y, data$Z, ifelse(data$Z == 0, data$w_observed_true, 1))
  } else {
    mean(data$Y[data$Z == 1]) - stats::weighted.mean(data$Y[idx0], w_obs0)
  }

  tau_ideal <- if (exists("estimate_att", mode = "function")) {
    estimate_att(data$Y, data$Z, ifelse(data$Z == 0, data$w_ideal_true, 1))
  } else {
    mean(data$Y[data$Z == 1]) - stats::weighted.mean(data$Y[idx0], w_ideal0)
  }

  ps_iqr <- as.numeric(stats::IQR(data$ps_observed_true, na.rm = TRUE))
  weight_cv <- ifelse(abs(mean_w) > 0, sqrt(var_observed) / abs(mean_w), NA_real_)

  denominator_small <- !is.finite(var_observed) || var_observed <= eps
  ratio_large <- is.infinite(original_ratio) ||
    (is.finite(original_ratio) && original_ratio >= ratio_warning_threshold)
  good_overlap_signature <- is.finite(ps_iqr) &&
    ps_iqr <= ps_iqr_warning_threshold &&
    is.finite(weight_cv) &&
    weight_cv <= weight_cv_warning_threshold

  flag <- denominator_small || ratio_large || good_overlap_signature

  status <- ifelse(
    !is.finite(var_observed) || var_observed <= variance_floor,
    "degenerate_or_near_degenerate_observed_weight_variance",
    ifelse(flag, "unstable_good_overlap_warning", "stable")
  )

  data.frame(
    observed_strength = observed_strength,
    hidden_strength = hidden_strength,
    n_total = nrow(data),
    n_control = sum(idx0),
    ps_observed_min = min(data$ps_observed_true, na.rm = TRUE),
    ps_observed_max = max(data$ps_observed_true, na.rm = TRUE),
    ps_observed_iqr = ps_iqr,
    weight_cv_control = weight_cv,
    var_observed_control = var_observed,
    var_ideal_control = var_ideal,
    epsilon_added_to_denominator = eps,
    original_variance_ratio = original_ratio,
    stabilized_variance_ratio = stabilized_ratio,
    original_R2 = original_R2,
    stabilized_R2 = stabilized_R2,
    tau_observed = tau_observed,
    tau_ideal = tau_ideal,
    tau_gap_observed_minus_ideal = tau_observed - tau_ideal,
    original_bound = original_bound,
    stabilized_bound = stabilized_bound,
    flag_good_overlap_instability = flag,
    diagnostic_status = status,
    stringsAsFactors = FALSE
  )
}

run_good_overlap_instability_grid <- function(
    observed_strength_grid = c(0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.4),
    replicates = 60,
    n = 1800,
    hidden_strength = 0.12,
    seed = 90210,
    epsilon_fraction = 1e-4,
    variance_floor = 1e-10
) {
  rows <- list()
  k <- 1
  for (a in observed_strength_grid) {
    for (r in seq_len(replicates)) {
      d <- make_good_overlap_instability_data(
        n = n,
        observed_strength = a,
        hidden_strength = hidden_strength,
        seed = seed + 1000 * match(a, observed_strength_grid) + r
      )
      row <- good_overlap_metric_row(
        data = d,
        observed_strength = a,
        hidden_strength = hidden_strength,
        epsilon_fraction = epsilon_fraction,
        variance_floor = variance_floor
      )
      row$replicate <- r
      rows[[k]] <- row
      k <- k + 1
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

summarize_good_overlap_grid <- function(replicate_metrics) {
  numeric_cols <- names(replicate_metrics)[
    vapply(replicate_metrics, function(x) is.numeric(x) && !is.logical(x), logical(1))
  ]
  keep_cols <- setdiff(numeric_cols, c("replicate", "n_total", "n_control", "hidden_strength"))

  pieces <- lapply(
    split(replicate_metrics, replicate_metrics$observed_strength),
    function(d) {
      values <- data.frame(
        observed_strength = unique(d$observed_strength),
        hidden_strength = unique(d$hidden_strength)[1],
        n_replicates = nrow(d),
        flag_rate = mean(d$flag_good_overlap_instability, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      for (nm in keep_cols) {
        x <- d[[nm]]
        finite_x <- x[is.finite(x)]
        values[[paste0(nm, "_median")]] <- ifelse(length(finite_x) > 0, stats::median(finite_x), ifelse(any(is.infinite(x)), Inf, NA_real_))
        values[[paste0(nm, "_q05")]] <- ifelse(length(finite_x) > 0, as.numeric(stats::quantile(finite_x, 0.05, names = FALSE)), ifelse(any(is.infinite(x)), Inf, NA_real_))
        values[[paste0(nm, "_q95")]] <- ifelse(length(finite_x) > 0, as.numeric(stats::quantile(finite_x, 0.95, names = FALSE)), ifelse(any(is.infinite(x)), Inf, NA_real_))
      }
      values
    }
  )
  out <- do.call(rbind, pieces)
  out <- out[order(out$observed_strength), , drop = FALSE]
  rownames(out) <- NULL
  out
}
