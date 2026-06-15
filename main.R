# ============================================================
# Main script
# Variance-based sensitivity analysis with optional integrated RGM extension
# ============================================================

# ------------------------------------------------------------
# Robust project-root-aware source loading
# ------------------------------------------------------------
# Version10.2 fix:
#   - works when the zip was extracted with or without an outer folder;
#   - works when source("path/to/main.R") is called from a parent folder;
#   - searches child folders such as final_project_version10/ if needed.

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
    candidates <- c(
      candidates,
      dirname(normalizePath(anchor_file[1], winslash = "/", mustWork = FALSE))
    )
  }

  stack_files <- v10_source_stack_files()
  if (length(stack_files) > 0) {
    candidates <- c(
      candidates,
      dirname(normalizePath(stack_files, winslash = "/", mustWork = FALSE))
    )
  }

  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) {
    candidates <- c(
      candidates,
      dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE))
    )
  }

  candidates <- c(candidates, normalizePath(getwd(), winslash = "/", mustWork = FALSE))
  candidates <- unique(candidates[nzchar(candidates)])

  # 1) Check each candidate and its parents.
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

  # 2) Check common child-folder layout after unzipping.
  for (start in candidates) {
    if (!dir.exists(start)) next
    child_dirs <- list.dirs(start, recursive = FALSE, full.names = TRUE)
    child_dirs <- c(child_dirs, list.dirs(start, recursive = TRUE, full.names = TRUE))
    child_dirs <- child_dirs[grepl("final_project|version10|version_10|v10", basename(child_dirs), ignore.case = TRUE)]
    for (child in unique(child_dirs)) {
      if (v10_is_project_root(child)) {
        return(normalizePath(child, winslash = "/", mustWork = FALSE))
      }
    }
  }

  stop(
    "Cannot determine project root. Please unzip the full archive and run source('main.R') ",
    "from the folder that contains main.R and functions/config.R."
  )
}

.project_file <- tryCatch(
  normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE),
  error = function(e) NA_character_
)

.project_dir <- v10_locate_project_root(.project_file)
.project_env <- environment()

setwd(.project_dir)

source_project_file <- function(relative_path) {
  path <- file.path(.project_dir, relative_path)
  if (!file.exists(path)) {
    stop(
      paste0(
        "Required source file missing: ", path,
        "\nProject root detected as: ", .project_dir,
        "\nPlease make sure the full version10.2 archive was extracted, not only main.R."
      )
    )
  }
  source(path, local = .project_env, chdir = TRUE)
}


source_project_file("functions/config.R")
source_project_file("functions/extension_utils.R")
source_project_file("functions/data_prep.R")
source_project_file("functions/ps_model.R")
source_project_file("functions/vbm_bounds.R")
source_project_file("functions/bootstrap_analysis.R")
source_project_file("functions/msm_analysis.R")
source_project_file("functions/benchmark_analysis.R")
source_project_file("functions/vbm_corr_benchmark.R")
source_project_file("functions/msm_qbal_benchmark.R")
source_project_file("functions/plotting.R")

# ------------------------------------------------------------
# Verify that all active project functions were loaded
# ------------------------------------------------------------

required_project_functions <- c(
  "prepare_data",
  "add_ps_model_columns_to_encoding_diagnostics",
  "fit_ps_model",
  "clip_ps",
  "compute_att_weight_function",
  "compute_weights",
  "estimate_att",
  "finite_sample_variance",
  "bias_scale_base",
  "bias_correlation_limit",
  "bias_scale",
  "bias_bound",
  "find_R2_bootstrap",
  "bootstrap_ci_from_stats_corr",
  "bootstrap_stats_summary",
  "generate_msm_bootstrap_curve",
  "msm_att_bounds",
  "msm_qbal_ratio_summary",
  "msm_qbal_att_bounds",
  "qbal_effective_gamma_from_benchmark",
  "compute_benchmark_R2",
  "compute_benchmark_corr",
  "compute_benchmark_Gamma",
  "make_vbm_corr_benchmark_row",
  "make_msm_qbal_benchmark_row",
  "generate_vbm_benchmark_bootstrap_curve",
  "generate_msm_benchmark_bootstrap_curve",
  "generate_benchmark_comparison_data",
  "run_covariate_benchmarks",
  "plot_bootstrap_curve",
  "normalize_benchmark_method_label",
  "plot_covariate_benchmark"
)

missing_project_functions <- required_project_functions[
  !vapply(
    required_project_functions,
    function(x) exists(x, mode = "function", envir = .project_env, inherits = FALSE),
    logical(1)
  )
]

if (length(missing_project_functions) > 0) {
  stop(
    paste0(
      "Missing required project function(s): ",
      paste(missing_project_functions, collapse = ", "),
      ". Run aborted to avoid unresolved project dependencies."
    )
  )
}

cat("Project self-check passed: all required project functions were loaded from ", .project_dir, ".\n", sep = "")

# ------------------------------------------------------------
# Prepare output folders
# ------------------------------------------------------------

dir.create(config$output_figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(config$output_tables_dir, recursive = TRUE, showWarnings = FALSE)

if (isTRUE(config$clear_output_before_run)) {
  old_tables <- list.files(
    config$output_tables_dir,
    pattern = "\\.csv$",
    full.names = TRUE
  )
  old_figures <- list.files(
    config$output_figures_dir,
    pattern = "\\.(png|pdf|svg)$",
    full.names = TRUE,
    ignore.case = TRUE
  )

  unlink(c(old_tables, old_figures), force = TRUE)

  cat(
    "Existing output CSV/PNG/PDF/SVG files removed before this run to keep results traceable.\n"
  )
}

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

if (!file.exists(config$data_path)) {
  stop(
    paste0(
      "Data file not found: ", config$data_path, "\n",
      "Please put nhanes.fish.rda into the data/ folder."
    )
  )
}

load(config$data_path)

if (!exists("nhanes.fish")) {
  stop("The loaded .rda file should contain an object named nhanes.fish.")
}

# ------------------------------------------------------------
# Prepare analysis data
# ------------------------------------------------------------

analysis <- prepare_data(nhanes.fish, config = config)

cat("Sample size:", nrow(analysis), "\n")
cat("Treated:", sum(analysis$Z == 1), "\n")
cat("Control:", sum(analysis$Z == 0), "\n")

ps_encoding_diagnostics <- attr(analysis, "encoding_diagnostics")
ps_encoding_diagnostics <- add_ps_model_columns_to_encoding_diagnostics(
  diagnostics = ps_encoding_diagnostics,
  full_formula = ps_formula,
  data = analysis,
  config = config
)

cat("PS covariate encoding:\n")
print(ps_encoding_diagnostics[, c(
  "variable",
  "final_encoding",
  "final_n_unique_after_drop",
  "model_n_columns",
  "model_encoding_status"
)])

# ------------------------------------------------------------
# Fit propensity score model and compute ATT weights
# ------------------------------------------------------------

set.seed(config$seed)

ps_fit <- fit_ps_model(
  formula = ps_formula,
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

# ------------------------------------------------------------
# ATT estimation
# ------------------------------------------------------------

tau_hat <- estimate_att(
  Y = analysis$Y,
  Z = analysis$Z,
  w = analysis$w
)

cat("ATT =", tau_hat, "\n")

# ------------------------------------------------------------
# Variance-based bootstrap analysis
# ------------------------------------------------------------
# The retained ATT-R2* curve uses percentile-bootstrap inference.

bootstrap_output <- find_R2_bootstrap(
  data = analysis,
  config = config
)

bootstrap_results <- bootstrap_output$results
R2_boot <- bootstrap_output$R2_star
vbm_bootstrap_stats_summary <- bootstrap_output$bootstrap_stats_summary

vbm_point_bounds_curve <- data.frame(
  R2 = seq(
    from = config$R2_min,
    to = config$R2_max,
    by = config$R2_coarse_step
  ),
  stringsAsFactors = FALSE
)

vbm_point_bounds_curve$bias_bound <- vapply(
  vbm_point_bounds_curve$R2,
  function(r) {
    bias_bound(
      w = analysis$w,
      Y = analysis$Y,
      Z = analysis$Z,
      R2 = r,
      config = config
    )
  },
  numeric(1)
)

vbm_point_bounds_curve$tau_hat <- tau_hat
vbm_point_bounds_curve$lower <- tau_hat - vbm_point_bounds_curve$bias_bound
vbm_point_bounds_curve$upper <- tau_hat + vbm_point_bounds_curve$bias_bound
vbm_point_bounds_curve$interval_type <- "closed_form_point_estimate_bounds"

cat("Bootstrap R²* =", R2_boot, "\n")

# ------------------------------------------------------------
# Covariate benchmark comparison
# ------------------------------------------------------------
# Figure 3 benchmark output contains four model rows per covariate:
#   MSM, MSM (Qbal), VBM, and VBM, w/ Corr.

benchmark_output <- run_covariate_benchmarks(
  analysis = analysis,
  tau_hat = tau_hat,
  full_formula = ps_formula,
  config = config
)

benchmark_df <- benchmark_output$benchmark_df
vbm_benchmark_bootstrap_curve <- benchmark_output$vbm_benchmark_bootstrap_curve
msm_benchmark_bootstrap_curve <- benchmark_output$msm_benchmark_bootstrap_curve

if (is.null(vbm_benchmark_bootstrap_curve)) {
  vbm_benchmark_bootstrap_curve <- data.frame()
}

if (is.null(msm_benchmark_bootstrap_curve)) {
  msm_benchmark_bootstrap_curve <- data.frame()
}

# ------------------------------------------------------------
# Save concise tables
# ------------------------------------------------------------
# Retained outputs cover bootstrap VBM inference and four-model benchmark data.

unweighted_att <- mean(analysis$Y[analysis$Z == 1]) -
  mean(analysis$Y[analysis$Z == 0])

paper_metrics <- c(
  "n_total",
  "n_treated",
  "n_control",
  "unweighted_ATT",
  "IPW_ATT",
  "VBM_R2_star_bootstrap",
  "ps_clip_lower",
  "ps_clip_upper",
  "weight_truncation"
)

paper_reported <- c(
  1107,
  234,
  873,
  2.37,
  2.14,
  0.52,
  NA_real_,
  NA_real_,
  NA_real_
)

current_version <- c(
  nrow(analysis),
  sum(analysis$Z == 1),
  sum(analysis$Z == 0),
  unweighted_att,
  tau_hat,
  R2_boot,
  config$ps_clip_lower,
  config$ps_clip_upper,
  config$weight_truncation
)

paper_notes <- c(
  "NHANES analysis sample size reported in the paper.",
  "Number of treated/high fish-shellfish consumers reported in the paper.",
  "Number of controls reported in the paper.",
  "Paper reports the unweighted ATT rounded to 2.37.",
  "Paper reports the IPW ATT rounded to 2.14.",
  "Paper reports the NHANES VBM threshold R2* around 0.52.",
  "Current implementation setting.",
  "Current implementation setting.",
  "Current implementation setting for main ATT weights."
)

final_results_vs_paper <- data.frame(
  metric = paper_metrics,
  paper_reported = paper_reported,
  current_version = current_version,
  difference_current_minus_paper = current_version - paper_reported,
  notes = paper_notes,
  stringsAsFactors = FALSE
)

benchmark_plot_data <- benchmark_df[
  ,
  intersect(
    c(
      "variable",
      "removed_terms",
      "method",
      "sensitivity_parameter",
      "sensitivity_value",
      "lower",
      "upper",
      "tau_hat",
      "interval_type",
      "inference_type",
      "bootstrap_B",
      "bootstrap_n_valid",
      "bootstrap_required_valid",
      "bootstrap_pool",
      "benchmark_status",
      "benchmark_group_used",

      "R2_minus_raw",
      "R2_directional",
      "R2_control",
      "R2_treated",
      "R2_all",
      "R2_control_directional",
      "R2_treated_directional",
      "R2_all_directional",
      "var_full",
      "var_reduced",

      "corr_raw",
      "corr_used",
      "corr_limit",
      "corr_group_used",

      "gamma_raw",
      "gamma_used",
      "gamma_group_used",
      "ratio_max",
      "ratio_q95",
      "ratio_q99",
      "ratio_qbal",
      "qbal_effective_gamma",
      "qbal_ratio_quantile",
      "n_ratio",

      "full_model_n_columns",
      "reduced_model_n_columns",
      "removed_model_columns"
    ),
    names(benchmark_df)
  )
]

vbm_bootstrap_plot_data <- bootstrap_results[
  ,
  intersect(
    c("stage", "R2", "lower", "upper", "n_valid", "B"),
    names(bootstrap_results)
  )
]

write.csv(
  final_results_vs_paper,
  file.path(config$output_tables_dir, "final_results_vs_paper.csv"),
  row.names = FALSE
)

write.csv(
  vbm_bootstrap_plot_data,
  file.path(config$output_tables_dir, "plot_vbm_bootstrap_curve.csv"),
  row.names = FALSE
)

write.csv(
  vbm_point_bounds_curve,
  file.path(config$output_tables_dir, "plot_vbm_point_bounds_curve.csv"),
  row.names = FALSE
)

write.csv(
  benchmark_plot_data,
  file.path(config$output_tables_dir, "plot_covariate_benchmark_vbm_msm.csv"),
  row.names = FALSE
)

write.csv(
  vbm_benchmark_bootstrap_curve,
  file.path(config$output_tables_dir, "plot_vbm_benchmark_bootstrap_curve.csv"),
  row.names = FALSE
)

write.csv(
  msm_benchmark_bootstrap_curve,
  file.path(config$output_tables_dir, "plot_msm_benchmark_bootstrap_curve.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Generate and save retained plots
# ------------------------------------------------------------

p1 <- plot_bootstrap_curve(
  results = bootstrap_results,
  R2_star = R2_boot
)

p2 <- plot_covariate_benchmark(
  df = benchmark_df,
  tau_hat = tau_hat,
  config = config
)

ggsave(
  file.path(config$output_figures_dir, "variance_based_bootstrap_curve.png"),
  p1,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(config$output_figures_dir, "covariate_benchmark_vbm_msm.png"),
  p2,
  width = 9.5,
  height = 5.5,
  dpi = 300
)


# ------------------------------------------------------------
# Main-program enhanced plots
# ------------------------------------------------------------
# This is the enhancedplot auto-launch for the base NHANES main program.
# It is separate from the four extension enhancedplot scripts.
if (isTRUE(config$run_main_enhanced_plots)) {
  main_enhanced_output_dir <- config$main_enhanced_output_dir
  if (is.null(main_enhanced_output_dir) || !nzchar(main_enhanced_output_dir)) {
    main_enhanced_output_dir <- file.path("output", "figures_enhanced")
  }

  v10_run_python_plot(
    project_dir = .project_dir,
    script_relative_path = "scripts/enhanced_plots.py",
    input_dir = config$output_tables_dir,
    output_dir = main_enhanced_output_dir,
    label = "Main-program enhanced plots"
  )
} else {
  cat("Main-program enhanced plots are disabled by config$run_main_enhanced_plots.\n")
}




# ------------------------------------------------------------
# Version10.2 root-level extension scripts
# ------------------------------------------------------------
# No extensions/ folder is required.  main.R directly sources the four
# standalone scripts below.  Each script can also be run independently from
# the project root:
#   source("extension_01_hidden_strength_vbm_msm.R")
#   source("extension_02_vbm_ps_misspecification.R")
#   source("extension_03_good_overlap_instability.R")
#   source("extension_04_rgm_vbm_msm.R")

run_v10_standalone_extension <- function(flag_name, script_name, extension_id, extension_label) {
  enabled <- isTRUE(config[[flag_name]])

  row <- data.frame(
    extension_id = extension_id,
    extension_label = extension_label,
    enabled_by_default = TRUE,
    enabled_this_run = enabled,
    entry_point = script_name,
    output_dir = file.path("output", "extension_results", extension_id),
    stringsAsFactors = FALSE
  )

  if (enabled) {
    cat("\n========== Running Version10.2 extension: ", extension_label, " ==========\n", sep = "")
    .v10_project_dir_from_runner <<- .project_dir
    source_project_file(script_name)
  } else {
    cat("Skipping Version10.2 extension ", extension_id, " because ", flag_name, " is FALSE.\n", sep = "")
  }

  row
}

if (isTRUE(config$run_version10_extensions)) {
  version10_extension_rows <- list(
    run_v10_standalone_extension(
      flag_name = "run_extension_01_hidden_strength_vbm_msm",
      script_name = "extension_01_hidden_strength_vbm_msm.R",
      extension_id = "01_hidden_strength_vbm_msm",
      extension_label = "Hidden-strength VBM/MSM comparison"
    ),
    run_v10_standalone_extension(
      flag_name = "run_extension_02_vbm_ps_misspecification",
      script_name = "extension_02_vbm_ps_misspecification.R",
      extension_id = "02_vbm_ps_misspecification",
      extension_label = "VBM under PS model misspecification"
    ),
    run_v10_standalone_extension(
      flag_name = "run_extension_03_good_overlap_instability",
      script_name = "extension_03_good_overlap_instability.R",
      extension_id = "03_good_overlap_instability",
      extension_label = "Good-overlap denominator instability"
    ),
    run_v10_standalone_extension(
      flag_name = "run_extension_04_rgm_vbm_msm",
      script_name = "extension_04_rgm_vbm_msm.R",
      extension_id = "04_rgm_vbm_msm",
      extension_label = "RGM model and VBM/MSM comparison"
    )
  )

  version10_index_dir <- file.path(.project_dir, "output", "extension_results")
  dir.create(version10_index_dir, recursive = TRUE, showWarnings = FALSE)
  version10_index <- do.call(rbind, version10_extension_rows)
  write.csv(
    version10_index,
    file.path(version10_index_dir, "extension_index.csv"),
    row.names = FALSE
  )
  cat("Version10.2 extension index saved to: ",
      file.path(version10_index_dir, "extension_index.csv"), "\n", sep = "")
} else {
  cat("Version10.2 extensions are disabled by config$run_version10_extensions.\n")
}


cat("\n========== Final results ==========" , "\n")
print(final_results_vs_paper)
cat("\nCleaned bootstrap and four-model benchmark outputs saved to output/figures and output/tables.\n")
