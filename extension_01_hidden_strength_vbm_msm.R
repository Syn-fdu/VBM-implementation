# ============================================================
# Version10 extension 01 main script
# Hidden-strength VBM/MSM comparison
# ============================================================

# This is root-level Version10.2 extension 01; it can be run by main.R or independently.
# Recommended run:
#   source("extension_01_hidden_strength_vbm_msm.R")
# or:
#   setwd(".../final_project_version7/synthetic_extension")
#   source("main_synthetic_vbm_msm.R")

# ------------------------------------------------------------
# Robust project-root detection
# ------------------------------------------------------------
v10_is_project_root <- function(path) {
  dir.exists(file.path(path, "functions")) &&
    file.exists(file.path(path, "functions", "config.R")) &&
    file.exists(file.path(path, "main.R"))
}

v10_source_stack_files <- function() {
  out <- character(0)
  for (i in seq_len(sys.nframe())) {
    candidate <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (!is.null(candidate) && length(candidate) > 0 && nzchar(candidate[1])) {
      out <- c(out, candidate[1])
    }
  }
  unique(out)
}

v10_locate_project_root <- function(anchor_file = NULL) {
  candidates <- character(0)

  if (!is.null(anchor_file) && length(anchor_file) > 0 && nzchar(anchor_file[1])) {
    candidates <- c(candidates, dirname(normalizePath(anchor_file[1], winslash = "/", mustWork = FALSE)))
  }

  stack_files <- v10_source_stack_files()
  if (length(stack_files) > 0) {
    candidates <- c(candidates, dirname(normalizePath(stack_files, winslash = "/", mustWork = FALSE)))
  }

  candidates <- c(candidates, normalizePath(getwd(), winslash = "/", mustWork = FALSE))
  candidates <- unique(candidates[nzchar(candidates)])

  for (start in candidates) {
    current <- start
    for (i in seq_len(12)) {
      if (v10_is_project_root(current)) {
        return(normalizePath(current, winslash = "/", mustWork = FALSE))
      }
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }

  for (start in candidates) {
    if (!dir.exists(start)) next
    child_dirs <- list.dirs(start, recursive = TRUE, full.names = TRUE)
    child_dirs <- child_dirs[grepl("final_project|version10|version_10|v10", basename(child_dirs), ignore.case = TRUE)]
    for (child in unique(child_dirs)) {
      if (v10_is_project_root(child)) {
        return(normalizePath(child, winslash = "/", mustWork = FALSE))
      }
    }
  }

  stop(
    "Cannot locate project root. Please unzip the full archive and run from the folder ",
    "containing main.R and functions/config.R."
  )
}

v10_detect_project_dir <- function() {
  if (exists(".v10_project_dir_from_runner")) {
    return(normalizePath(.v10_project_dir_from_runner, winslash = "/", mustWork = FALSE))
  }
  source_file <- tryCatch(
    normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE),
    error = function(e) NA_character_
  )
  v10_locate_project_root(source_file)
}

.project_dir <- v10_detect_project_dir()

.project_env <- environment()
setwd(.project_dir)

source_project_file <- function(relative_path) {
  path <- file.path(.project_dir, relative_path)
  if (!file.exists(path)) {
    stop(paste0("Required source file missing: ", path))
  }
  source(path, local = .project_env, chdir = FALSE)
}

source_project_file("functions/config.R")
source_project_file("functions/ps_model.R")
source_project_file("functions/vbm_bounds.R")
source_project_file("functions/bootstrap_analysis.R")
source_project_file("functions/msm_analysis.R")
source_project_file("functions/benchmark_analysis.R")
source_project_file("functions/vbm_corr_benchmark.R")
source_project_file("functions/msm_qbal_benchmark.R")
source_project_file("functions/plotting.R")
source_project_file("functions/synthetic_data_generation.R")

cat("Synthetic project root: ", .project_dir, "\n", sep = "")

# ------------------------------------------------------------
# Synthetic-specific config
# ------------------------------------------------------------
config$seed <- 2027

config$output_figures_dir <- "output/extension_results/01_hidden_strength_vbm_msm/figures"
config$output_tables_dir <- "output/extension_results/01_hidden_strength_vbm_msm/tables"
dir.create(config$output_figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(config$output_tables_dir, recursive = TRUE, showWarnings = FALSE)

# Faster synthetic bootstrap settings.
config$R2_min <- 0
config$R2_max <- 0.90
config$R2_coarse_step <- 0.10
config$R2_fine_step <- 0.02
config$R2_fine_window <- 0.08
config$coarse_B <- 80
config$fine_B <- 120
config$vbm_benchmark_bootstrap_B <- 120
config$msm_benchmark_bootstrap_B <- 120
config$msm_bootstrap_B <- 120
config$diagnostics <- TRUE
config$min_bootstrap_valid_fraction <- 0.65
config$min_bootstrap_valid <- 25
config$weight_truncation <- 0.996
config$benchmark_weight_truncation <- NA_real_

# Keep synthetic benchmark table focused on the two models requested for
# Version10 extension 01: VBM benchmark and MSM benchmark across hidden strengths.
config$figure3_include_vbm_corr <- FALSE
config$figure3_include_qbal <- FALSE

observed_ps_formula <- Z ~ gender + age + income +
  income.missing + race +
  education + smoking.ever +
  smoking.now

base_benchmark_groups <- list(
  gender = c("gender"),
  age = c("age"),
  income = c("income"),
  income_missing = c("income.missing"),
  education = c("education"),
  cig_smoked = c("smoking.now"),
  smoking_history = c("smoking.ever"),
  race = c("race")
)

config$benchmark_groups <- base_benchmark_groups
config$benchmark_plot_labels <- c(
  gender = "Gender",
  age = "Age",
  income = "Income",
  income_missing = "Income (Missing)",
  education = "Education",
  cig_smoked = "Cig. Smoked",
  smoking_history = "Smoking history",
  race = "Race",
  location_region = "Location (oracle)",
  hidden_confounder = "Hidden confounder (oracle)"
)

# Remove stale synthetic outputs only.  The main NHANES output folders are not
# touched, which keeps source("main.R") isolated from the extension workflow.
clear_synthetic_outputs <- function(config, project_dir) {
  dirs <- c(
    file.path(project_dir, config$output_figures_dir),
    file.path(project_dir, config$output_tables_dir),
    file.path(project_dir, "output/extension_results/01_hidden_strength_vbm_msm/figures_enhanced")
  )

  for (d in dirs) {
    if (!dir.exists(d)) next
    files <- list.files(
      d,
      pattern = "\\.(csv|png|svg|pdf|jpg|jpeg|webp)$",
      full.names = TRUE,
      ignore.case = TRUE
    )
    if (length(files) > 0) unlink(files)
  }
}

synthetic_param <- function(analysis, name, default = NA) {
  params <- attr(analysis, "synthetic_params")
  if (!is.null(params) && !is.null(params[[name]])) {
    return(params[[name]])
  }

  val <- attr(analysis, name)
  if (!is.null(val)) {
    return(val)
  }

  default
}

clear_synthetic_outputs(config = config, project_dir = .project_dir)

# bootstrap_analysis.R uses the global object ps_formula.
ps_formula <- observed_ps_formula

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
estimate_ipw_att_from_formula <- function(analysis, formula, config) {
  fit <- fit_ps_model(
    formula = formula,
    data = analysis,
    config = config
  )

  ps <- clip_ps(
    ps = fit$fitted.values,
    lower = config$ps_clip_lower,
    upper = config$ps_clip_upper
  )

  w <- compute_weights(
    ps = ps,
    Z = analysis$Z,
    config = config
  )

  estimate_att(
    Y = analysis$Y,
    Z = analysis$Z,
    w = w
  )
}

make_oracle_benchmark_rows <- function(
    analysis,
    tau_hat,
    oracle_formula,
    omitted_variable,
    omitted_terms,
    config
) {

  r2_info <- compute_benchmark_R2(
    full_formula = oracle_formula,
    benchmark_name = omitted_variable,
    remove_terms = omitted_terms,
    data = analysis,
    config = config
  )

  r2 <- r2_info$R2_used

  vbm_curve <- generate_vbm_benchmark_bootstrap_curve(
    analysis = analysis,
    R2_values = c(r2),
    config = config
  )

  vbm_row <- select_nearest_R2_row(
    vbm_bootstrap_curve = vbm_curve,
    R2 = r2
  )

  bb <- bias_bound(
    w = analysis$w,
    Y = analysis$Y,
    Z = analysis$Z,
    R2 = r2,
    config = config
  )

  det_lower <- tau_hat - bb
  det_upper <- tau_hat + bb

  if (!is.null(vbm_row) && is.finite(vbm_row$lower) && is.finite(vbm_row$upper)) {
    vbm_lower <- vbm_row$lower
    vbm_upper <- vbm_row$upper
    interval_type <- "bootstrap_percentile_fixed_R2"
    inference_type <- "bootstrap"
    boot_lower <- vbm_row$lower
    boot_upper <- vbm_row$upper
    boot_B <- vbm_row$B_requested
    boot_n <- vbm_row$n_valid
    boot_req <- vbm_row$required_valid
    boot_pool <- vbm_row$bootstrap_pool
  } else {
    vbm_lower <- det_lower
    vbm_upper <- det_upper
    interval_type <- "closed_form_fallback_no_valid_benchmark_bootstrap"
    inference_type <- "closed_form_fallback"
    boot_lower <- NA_real_
    boot_upper <- NA_real_
    boot_B <- NA_real_
    boot_n <- NA_real_
    boot_req <- NA_real_
    boot_pool <- NA_character_
  }

  out_vbm <- make_benchmark_row(
    variable = omitted_variable,
    removed_terms = format_removed_terms(omitted_terms),
    method = "VBM benchmark (bootstrap)",
    sensitivity_parameter = "R2",
    sensitivity_value = r2,
    lower = vbm_lower,
    upper = vbm_upper,
    interval_type = interval_type,
    inference_type = inference_type,
    tau_hat = tau_hat,
    deterministic_lower = det_lower,
    deterministic_upper = det_upper,
    bootstrap_lower = boot_lower,
    bootstrap_upper = boot_upper,
    bootstrap_B = boot_B,
    bootstrap_n_valid = boot_n,
    bootstrap_required_valid = boot_req,
    bootstrap_pool = boot_pool,
    benchmark_status = paste(
      r2_info$benchmark_status,
      "oracle_omitted_variable_benchmark",
      sep = " | "
    ),
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

  gamma_info <- compute_benchmark_Gamma(
    full_formula = oracle_formula,
    benchmark_name = omitted_variable,
    remove_terms = omitted_terms,
    data = analysis,
    config = config
  )

  gamma <- gamma_info$gamma_used

  msm_curve <- generate_msm_benchmark_bootstrap_curve(
    analysis = analysis,
    gamma_values = c(gamma),
    config = config
  )

  msm_row <- select_nearest_gamma_row(
    msm_bootstrap_curve = msm_curve,
    gamma = gamma
  )

  det_msm <- msm_att_bounds(
    analysis = analysis,
    Gamma = gamma,
    config = config
  )

  if (!is.null(msm_row) && is.finite(msm_row$lower) && is.finite(msm_row$upper)) {
    msm_lower <- msm_row$lower
    msm_upper <- msm_row$upper
    msm_interval_type <- "bootstrap_percentile_fixed_Gamma"
    msm_inference_type <- "bootstrap"
    msm_boot_lower <- msm_row$lower
    msm_boot_upper <- msm_row$upper
    msm_boot_B <- msm_row$B_requested
    msm_boot_n <- msm_row$n_valid
    msm_boot_req <- msm_row$required_valid
    msm_boot_pool <- msm_row$bootstrap_pool
  } else {
    msm_lower <- det_msm$lower
    msm_upper <- det_msm$upper
    msm_interval_type <- "closed_form_fallback_no_valid_benchmark_bootstrap"
    msm_inference_type <- "closed_form_fallback"
    msm_boot_lower <- NA_real_
    msm_boot_upper <- NA_real_
    msm_boot_B <- NA_real_
    msm_boot_n <- NA_real_
    msm_boot_req <- NA_real_
    msm_boot_pool <- NA_character_
  }

  out_msm <- make_benchmark_row(
    variable = omitted_variable,
    removed_terms = format_removed_terms(omitted_terms),
    method = "MSM benchmark (bootstrap)",
    sensitivity_parameter = "Gamma",
    sensitivity_value = gamma,
    lower = msm_lower,
    upper = msm_upper,
    interval_type = msm_interval_type,
    inference_type = msm_inference_type,
    tau_hat = tau_hat,
    deterministic_lower = det_msm$lower,
    deterministic_upper = det_msm$upper,
    bootstrap_lower = msm_boot_lower,
    bootstrap_upper = msm_boot_upper,
    bootstrap_B = msm_boot_B,
    bootstrap_n_valid = msm_boot_n,
    bootstrap_required_valid = msm_boot_req,
    bootstrap_pool = msm_boot_pool,
    benchmark_status = paste(
      gamma_info$status,
      "oracle_omitted_variable_benchmark",
      sep = " | "
    ),
    gamma_raw = gamma_info$gamma_raw,
    gamma_used = gamma_info$gamma_used,
    gamma_group_used = gamma_info$gamma_group_used,
    ratio_max = gamma_info$ratio_max,
    ratio_q95 = gamma_info$ratio_q95,
    ratio_q99 = gamma_info$ratio_q99,
    ratio_qbal = gamma_info$ratio_qbal,
    n_ratio = gamma_info$n_ratio
  )

  out <- rbind(out_vbm, out_msm)
  out$benchmark_role <- "oracle_omitted_confounder"
  out
}


save_synthetic_standard_plots <- function(
    bootstrap_results,
    benchmark_df,
    tau_hat,
    R2_star,
    output_prefix,
    dataset_label,
    config
) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cat("Skipping standard synthetic plots because ggplot2 is not available.\n")
    return(invisible(FALSE))
  }

  dataset_label <- as.character(dataset_label)[1]

  dir.create(config$output_figures_dir, recursive = TRUE, showWarnings = FALSE)

  curve_png <- file.path(
    config$output_figures_dir,
    paste0(output_prefix, "_vbm_bootstrap_curve.png")
  )

  benchmark_png <- file.path(
    config$output_figures_dir,
    paste0(output_prefix, "_benchmark_vbm_msm.png")
  )

  p_curve <- plot_bootstrap_curve(
    results = bootstrap_results,
    R2_star = R2_star
  ) +
    ggplot2::labs(
      title = paste0("Synthetic VBM bootstrap curve: ", dataset_label),
      subtitle = "Vertical dotted line marks the searched R2* threshold.",
      x = expression(R^2),
      y = "Bootstrap ATT interval"
    )

  tryCatch(
    {
      ggplot2::ggsave(
        filename = curve_png,
        plot = p_curve,
        width = 10.5,
        height = 6.2,
        dpi = 300
      )
      cat("Saved VBM bootstrap curve plot: ", curve_png, "\n", sep = "")
    },
    error = function(e) {
      cat("Failed to save VBM bootstrap curve plot for ", output_prefix, ": ",
          conditionMessage(e), "\n", sep = "")
    }
  )

  p_benchmark <- plot_covariate_benchmark(
    df = benchmark_df,
    tau_hat = tau_hat,
    config = config
  ) +
    ggplot2::labs(
      title = paste0("Synthetic benchmark comparison: ", dataset_label),
      subtitle = "Known-covariate deletion benchmarks plus the oracle omitted-confounder benchmark.",
      x = "Benchmarked covariate",
      y = "ATT interval",
      color = "Sensitivity model"
    )

  tryCatch(
    {
      ggplot2::ggsave(
        filename = benchmark_png,
        plot = p_benchmark,
        width = 12.5,
        height = 6.8,
        dpi = 300
      )
      cat("Saved benchmark comparison plot: ", benchmark_png, "\n", sep = "")
    },
    error = function(e) {
      cat("Failed to save benchmark comparison plot for ", output_prefix, ": ",
          conditionMessage(e), "\n", sep = "")
    }
  )

  invisible(TRUE)
}

run_one_synthetic_dataset <- function(
    analysis,
    oracle_formula,
    omitted_variable,
    omitted_terms,
    output_prefix,
    config
) {

  ps_formula <<- observed_ps_formula

  config$benchmark_groups <- base_benchmark_groups

  ps_fit <- fit_ps_model(
    formula = observed_ps_formula,
    data = analysis,
    config = config
  )

  analysis$ps <- clip_ps(
    ps = ps_fit$fitted.values,
    lower = config$ps_clip_lower,
    upper = config$ps_clip_upper
  )

  analysis$w <- compute_weights(
    ps = analysis$ps,
    Z = analysis$Z,
    config = config
  )

  tau_hat <- estimate_att(
    Y = analysis$Y,
    Z = analysis$Z,
    w = analysis$w
  )

  unweighted_att <- mean(analysis$Y[analysis$Z == 1]) -
    mean(analysis$Y[analysis$Z == 0])

  oracle_att <- estimate_ipw_att_from_formula(
    analysis = analysis,
    formula = oracle_formula,
    config = config
  )

  true_att <- mean(analysis$true_unit_ATT[analysis$Z == 1])

  cat("\n========== Running ", unique(analysis$dataset_label), " ==========\n", sep = "")
  cat("n = ", nrow(analysis), "; treated = ", sum(analysis$Z == 1), "\n", sep = "")
  cat("True ATT among treated = ", round(true_att, 4), "\n", sep = "")
  cat("Observed-only IPW ATT = ", round(tau_hat, 4), "\n", sep = "")
  cat("Oracle omitted-variable IPW ATT = ", round(oracle_att, 4), "\n", sep = "")

  bootstrap_output <- find_R2_bootstrap(
    data = analysis,
    config = config
  )

  bootstrap_results <- bootstrap_output$results
  R2_star <- bootstrap_output$R2_star

  bootstrap_results$dataset_id <- unique(analysis$dataset_id)
  bootstrap_results$dataset_label <- unique(analysis$dataset_label)

  benchmark_output <- run_covariate_benchmarks(
    analysis = analysis,
    tau_hat = tau_hat,
    full_formula = observed_ps_formula,
    config = config
  )

  benchmark_df <- benchmark_output$benchmark_df
  benchmark_df$benchmark_role <- "observed_covariate"

  oracle_rows <- make_oracle_benchmark_rows(
    analysis = analysis,
    tau_hat = tau_hat,
    oracle_formula = oracle_formula,
    omitted_variable = omitted_variable,
    omitted_terms = omitted_terms,
    config = config
  )

  benchmark_df <- rbind(benchmark_df, oracle_rows)

  benchmark_df$dataset_id <- unique(analysis$dataset_id)
  benchmark_df$dataset_label <- unique(analysis$dataset_label)

  vbm_oracle_r2 <- benchmark_df$sensitivity_value[
    benchmark_df$variable == omitted_variable &
      benchmark_df$method == "VBM benchmark (bootstrap)"
  ][1]

  msm_oracle_gamma <- benchmark_df$sensitivity_value[
    benchmark_df$variable == omitted_variable &
      benchmark_df$method == "MSM benchmark (bootstrap)"
  ][1]

  known_vbm <- benchmark_df[
    benchmark_df$method == "VBM benchmark (bootstrap)" &
      benchmark_df$benchmark_role == "observed_covariate",
    ,
    drop = FALSE
  ]

  max_known_vbm_r2 <- suppressWarnings(max(known_vbm$sensitivity_value, na.rm = TRUE))
  if (!is.finite(max_known_vbm_r2)) {
    max_known_vbm_r2 <- NA_real_
  }

  # A scenario is flagged when the actually omitted benchmark is at least as
  # large as the searched R2* needed to cross the null. This directly
  # implements the paper's benchmark interpretation in the synthetic setting.
  overturn_success <- is.finite(R2_star) &&
    is.finite(vbm_oracle_r2) &&
    vbm_oracle_r2 >= R2_star

  dataset_summary <- data.frame(
    dataset_id = unique(analysis$dataset_id),
    dataset_label = unique(analysis$dataset_label),
    n_total = nrow(analysis),
    n_treated = sum(analysis$Z == 1),
    n_control = sum(analysis$Z == 0),
    treated_fraction = mean(analysis$Z == 1),
    omitted_variable = omitted_variable,
    scenario_order = synthetic_param(analysis, "scenario_order", NA_integer_),
    scenario_class = synthetic_param(analysis, "scenario_class", NA_character_),
    omitted_confounding_strength = synthetic_param(analysis, "omitted_confounding_strength", NA_real_),
    design_true_att = synthetic_param(analysis, "true_att", NA_real_),
    design_sigma_y = synthetic_param(analysis, "sigma_y", NA_real_),
    design_location_trt_strength = synthetic_param(analysis, "location_trt_strength", NA_real_),
    design_location_y_strength = synthetic_param(analysis, "location_y_strength", NA_real_),
    design_hidden_trt_strength = synthetic_param(analysis, "hidden_trt_strength", NA_real_),
    design_hidden_y_strength = synthetic_param(analysis, "hidden_y_strength", NA_real_),
    true_ATT_among_treated = true_att,
    unweighted_ATT = unweighted_att,
    observed_only_IPW_ATT = tau_hat,
    oracle_omitted_IPW_ATT = oracle_att,
    observed_minus_oracle_IPW_ATT = tau_hat - oracle_att,
    observed_minus_true_ATT = tau_hat - true_att,
    oracle_minus_true_ATT = oracle_att - true_att,
    VBM_R2_star_observed_only = R2_star,
    oracle_omitted_VBM_benchmark_R2 = vbm_oracle_r2,
    max_observed_covariate_VBM_benchmark_R2 = max_known_vbm_r2,
    oracle_omitted_MSM_benchmark_Gamma = msm_oracle_gamma,
    benchmark_says_omitted_can_overturn = overturn_success,
    bootstrap_rows = nrow(bootstrap_results),
    stringsAsFactors = FALSE
  )

  # Backward-compatible aliases used by earlier A/B scripts.
  dataset_summary$oracle_location_IPW_ATT <- ifelse(
    omitted_variable == "location_region",
    oracle_att,
    NA_real_
  )
  dataset_summary$oracle_hidden_IPW_ATT <- ifelse(
    omitted_variable == "hidden_confounder",
    oracle_att,
    NA_real_
  )

  write.csv(
    analysis,
    file.path(config$output_tables_dir, paste0(output_prefix, "_analysis_data.csv")),
    row.names = FALSE
  )

  write.csv(
    bootstrap_results,
    file.path(config$output_tables_dir, paste0(output_prefix, "_vbm_bootstrap_curve.csv")),
    row.names = FALSE
  )

  write.csv(
    benchmark_df,
    file.path(config$output_tables_dir, paste0(output_prefix, "_benchmark_vbm_msm.csv")),
    row.names = FALSE
  )

  save_synthetic_standard_plots(
    bootstrap_results = bootstrap_results,
    benchmark_df = benchmark_df,
    tau_hat = tau_hat,
    R2_star = R2_star,
    output_prefix = output_prefix,
    dataset_label = unique(analysis$dataset_label),
    config = config
  )

  list(
    analysis = analysis,
    tau_hat = tau_hat,
    oracle_att = oracle_att,
    true_att = true_att,
    R2_star = R2_star,
    bootstrap_results = bootstrap_results,
    benchmark_df = benchmark_df,
    dataset_summary = dataset_summary
  )
}

# ------------------------------------------------------------
# Generate and run the current synthetic datasets
# ------------------------------------------------------------
synthetic_datasets <- generate_all_synthetic_datasets(n = 540)

runs <- list()

for (dataset_id in names(synthetic_datasets)) {
  analysis_i <- synthetic_datasets[[dataset_id]]
  omitted_variable_i <- synthetic_param(
    analysis_i,
    "omitted_variable",
    "hidden_confounder"
  )
  omitted_terms_i <- synthetic_param(
    analysis_i,
    "omitted_terms",
    omitted_variable_i
  )

  oracle_formula_i <- update(
    observed_ps_formula,
    paste(". ~ . +", paste(omitted_terms_i, collapse = " + "))
  )

  runs[[dataset_id]] <- run_one_synthetic_dataset(
    analysis = analysis_i,
    oracle_formula = oracle_formula_i,
    omitted_variable = omitted_variable_i,
    omitted_terms = omitted_terms_i,
    output_prefix = dataset_id,
    config = config
  )
}

all_summaries <- do.call(
  rbind,
  lapply(runs, function(x) x$dataset_summary)
)

all_benchmark <- do.call(
  rbind,
  lapply(runs, function(x) x$benchmark_df)
)

all_bootstrap <- do.call(
  rbind,
  lapply(runs, function(x) x$bootstrap_results)
)

model_label_from_method <- function(method) {
  ifelse(
    grepl("^VBM", method),
    "VBM",
    ifelse(grepl("^MSM", method), "MSM", as.character(method))
  )
}

make_detection_power_rows <- function(runs) {

  rows <- lapply(runs, function(x) {

    s <- x$dataset_summary[1, , drop = FALSE]
    b <- x$benchmark_df[
      x$benchmark_df$benchmark_role == "oracle_omitted_confounder" &
        x$benchmark_df$method %in% c("VBM benchmark (bootstrap)", "MSM benchmark (bootstrap)"),
      ,
      drop = FALSE
    ]

    if (nrow(b) == 0) {
      return(NULL)
    }

    out <- data.frame()

    for (i in seq_len(nrow(b))) {
      lower <- as.numeric(b$lower[i])
      upper <- as.numeric(b$upper[i])
      tau_hat <- as.numeric(s$observed_only_IPW_ATT[1])

      if (is.finite(tau_hat) && tau_hat < 0) {
        detection_margin <- upper
        null_reached <- is.finite(upper) && upper >= 0
      } else {
        detection_margin <- -lower
        null_reached <- is.finite(lower) && lower <= 0
      }

      null_in_interval <- is.finite(lower) && is.finite(upper) && lower <= 0 && upper >= 0

      out <- rbind(
        out,
        data.frame(
          dataset_id = as.character(s$dataset_id[1]),
          dataset_label = as.character(s$dataset_label[1]),
          scenario_order = as.integer(s$scenario_order[1]),
          scenario_class = as.character(s$scenario_class[1]),
          omitted_variable = as.character(s$omitted_variable[1]),
          omitted_confounding_strength = as.numeric(s$omitted_confounding_strength[1]),
          model = model_label_from_method(as.character(b$method[i])),
          method = as.character(b$method[i]),
          sensitivity_parameter = as.character(b$sensitivity_parameter[i]),
          sensitivity_value = as.numeric(b$sensitivity_value[i]),
          interval_lower = lower,
          interval_upper = upper,
          interval_width = upper - lower,
          tau_hat = tau_hat,
          true_ATT = as.numeric(s$true_ATT_among_treated[1]),
          oracle_omitted_IPW_ATT = as.numeric(s$oracle_omitted_IPW_ATT[1]),
          null_in_interval = null_in_interval,
          null_reached_or_crossed = null_reached,
          detection_margin = detection_margin,
          detection_status = ifelse(
            null_reached,
            "overturn detected at oracle benchmark",
            "not enough confounding at oracle benchmark"
          ),
          stringsAsFactors = FALSE
        )
      )
    }

    out
  })

  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) {
    return(data.frame())
  }

  out <- do.call(rbind, rows)
  out[order(out$scenario_order, out$model), , drop = FALSE]
}

model_detection_power <- make_detection_power_rows(runs)

write.csv(
  all_summaries,
  file.path(config$output_tables_dir, "synthetic_dataset_summary.csv"),
  row.names = FALSE
)

write.csv(
  all_benchmark,
  file.path(config$output_tables_dir, "all_synthetic_benchmark_vbm_msm.csv"),
  row.names = FALSE
)

write.csv(
  all_bootstrap,
  file.path(config$output_tables_dir, "all_synthetic_vbm_bootstrap_curves.csv"),
  row.names = FALSE
)

write.csv(
  model_detection_power,
  file.path(config$output_tables_dir, "synthetic_model_detection_power.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Basic ATT diagnostic plot
# ------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {

  att_plot_df <- do.call(
    rbind,
    lapply(runs, function(x) {
      data.frame(
        dataset_id = unique(x$analysis$dataset_id),
        dataset_label = unique(x$analysis$dataset_label),
        quantity = c("Observed-only IPW", "Oracle omitted-variable IPW", "True ATT"),
        att = c(x$tau_hat, x$oracle_att, x$true_att),
        stringsAsFactors = FALSE
      )
    })
  )

  att_plot_df$quantity <- factor(
    att_plot_df$quantity,
    levels = c("Observed-only IPW", "Oracle omitted-variable IPW", "True ATT")
  )

  p_att <- ggplot2::ggplot(
    att_plot_df,
    ggplot2::aes(x = dataset_label, y = att, shape = quantity)
  ) +
    ggplot2::geom_point(size = 3.4) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)
    ) +
    ggplot2::labs(
      title = "Synthetic extension: true and estimated ATT",
      subtitle = "Oracle estimate includes the deliberately omitted variable in the PS model.",
      x = "Synthetic dataset",
      y = "ATT on log2 mercury scale",
      shape = "Quantity"
    )

  ggplot2::ggsave(
    file.path(config$output_figures_dir, "synthetic_att_summary.png"),
    p_att,
    width = 11,
    height = 6,
    dpi = 300
  )
}


run_enhanced_hidden_strength_plots <- function(config, project_dir) {

  script_path <- file.path(project_dir, "scripts", "enhancedplot_hidden_strength.py")

  if (!file.exists(script_path)) {
    cat("Skipping enhanced hidden-strength plots because script is missing: ",
        script_path, "\n", sep = "")
    return(invisible(FALSE))
  }

  python_bin <- Sys.which("python")
  if (!nzchar(python_bin)) {
    python_bin <- Sys.which("python3")
  }

  if (!nzchar(python_bin)) {
    cat("Skipping enhanced hidden-strength plots because neither python nor python3 was found.\n")
    cat("You can still run the script manually after installing Python packages pandas and matplotlib.\n")
    return(invisible(FALSE))
  }

  input_dir <- normalizePath(
    file.path(project_dir, config$output_tables_dir),
    winslash = "/",
    mustWork = FALSE
  )

  output_dir <- normalizePath(
    file.path(project_dir, "output/extension_results/01_hidden_strength_vbm_msm/figures_enhanced"),
    winslash = "/",
    mustWork = FALSE
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cat("\nGenerating enhanced hidden-strength plots with Python...\n")
  cat("Python: ", python_bin, "\n", sep = "")
  cat("Script: ", script_path, "\n", sep = "")

  result <- tryCatch(
    {
      system2(
        command = python_bin,
        args = c(script_path, "--input", input_dir, "--output", output_dir),
        stdout = TRUE,
        stderr = TRUE
      )
    },
    warning = function(w) {
      cat("Enhanced hidden-strength plot warning: ", conditionMessage(w), "\n", sep = "")
      return(NULL)
    },
    error = function(e) {
      cat("Enhanced hidden-strength plot error: ", conditionMessage(e), "\n", sep = "")
      return(NULL)
    }
  )

  if (!is.null(result) && length(result) > 0) {
    cat(paste(result, collapse = "\n"), "\n", sep = "")
  }

  status <- attr(result, "status")
  if (!is.null(status) && status != 0) {
    cat("Enhanced hidden-strength plot script exited with status ", status, ".\n", sep = "")
    return(invisible(FALSE))
  }

  cat("Enhanced hidden-strength plots saved under: ", output_dir, "\n", sep = "")
  invisible(TRUE)
}

run_enhanced_hidden_strength_plots(config = config, project_dir = .project_dir)

cat("\n========== Synthetic extension summary ==========\n")
print(all_summaries)
cat("\nOutputs saved under: ", config$output_tables_dir, " and ", config$output_figures_dir, "\n", sep = "")
cat("Standard per-dataset R plots are saved under output/extension_results/01_hidden_strength_vbm_msm/figures/.\n")
cat("Enhanced Figure 1-style and benchmark plots are saved under output/extension_results/01_hidden_strength_vbm_msm/figures_enhanced/.\n")
cat("If Python was unavailable during the R run, generate enhanced plots manually with:\n")
cat("  python scripts/enhancedplot_hidden_strength.py --input output/extension_results/01_hidden_strength_vbm_msm/tables --output output/extension_results/01_hidden_strength_vbm_msm/figures_enhanced\n")
