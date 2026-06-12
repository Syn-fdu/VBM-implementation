# ============================================================
# Version10 Extension 04:
# RGM model exploration and comparison with MSM/VBM
# ============================================================
# Standalone run:
#   source("extension_04_rgm_vbm_msm.R")
#
# RGM sensitivity parameter:
#   T = 1/2 * sum_i |q_i - p_i|,
# where p is the observed normalized control-weight distribution and q is an
# admissible ideal control-weight distribution. T is total variation distance.
#
# This extension runs the sharp finite-sample RGM search, adds RGM benchmark
# rows to the VBM/MSM benchmark table, and saves a comparison guide.
# ============================================================

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

source(file.path(.project_dir, "functions", "extension_utils.R"), local = .project_env)
v10_source_files(
  project_dir = .project_dir,
  files = c(
    "functions/config.R",
    "functions/data_prep.R",
    "functions/ps_model.R",
    "functions/vbm_bounds.R",
    "functions/bootstrap_analysis.R",
    "functions/msm_analysis.R",
    "functions/benchmark_analysis.R",
    "functions/vbm_corr_benchmark.R",
    "functions/msm_qbal_benchmark.R",
    "functions/plotting.R",
    "rgm_model/rgm_conservative.R",
    "rgm_model/rgm_sharp.R",
    "rgm_model/rgm_benchmark.R"
  ),
  env = .project_env
)

dirs <- v10_prepare_extension_dirs(.project_dir, "04_rgm_vbm_msm")
v10_clear_output_files(dirs$tables, dirs$figures, dirs$enhanced)

config$seed <- 2024
config$output_figures_dir <- file.path("output", "extension_results", "04_rgm_vbm_msm", "figures")
config$output_tables_dir <- file.path("output", "extension_results", "04_rgm_vbm_msm", "tables")

# Runtime-aware settings for automatic extension execution.
config$coarse_B <- 120
config$fine_B <- 180
config$vbm_benchmark_bootstrap_B <- 120
config$msm_benchmark_bootstrap_B <- 120
config$msm_bootstrap_B <- 120
config$rgm_sharp_coarse_B <- 120
config$rgm_sharp_fine_B <- 180
config$rgm_benchmark_bootstrap_B <- 120
config$min_bootstrap_valid_fraction <- 0.60
config$min_bootstrap_valid <- 25
config$diagnostics <- TRUE
config$figure3_include_vbm_corr <- FALSE
config$figure3_include_qbal <- FALSE

if (!file.exists(config$data_path)) {
  stop(paste0("Data file not found: ", config$data_path))
}
load(config$data_path)
if (!exists("nhanes.fish")) stop("The loaded .rda file should contain nhanes.fish.")

analysis <- prepare_data(nhanes.fish, config = config)

ps_fit <- fit_ps_model(formula = ps_formula, data = analysis, config = config)
analysis$ps <- clip_ps(ps_fit$fitted.values, config$ps_clip_lower, config$ps_clip_upper)
analysis$w <- compute_weights(ps = analysis$ps, Z = analysis$Z, config = config)

tau_hat <- estimate_att(analysis$Y, analysis$Z, analysis$w)

cat("Extension 04 NHANES analysis sample size: ", nrow(analysis), "\n", sep = "")
cat("Extension 04 ATT = ", tau_hat, "\n", sep = "")

vbm_search <- find_R2_bootstrap(data = analysis, config = config)
vbm_R2_star <- vbm_search$R2_star
vbm_bootstrap_curve <- vbm_search$results

benchmark_output <- run_covariate_benchmarks(
  analysis = analysis,
  tau_hat = tau_hat,
  full_formula = ps_formula,
  config = config
)

benchmark_df <- benchmark_output$benchmark_df
vbm_benchmark_bootstrap_curve <- benchmark_output$vbm_benchmark_bootstrap_curve
msm_benchmark_bootstrap_curve <- benchmark_output$msm_benchmark_bootstrap_curve

if (is.null(vbm_benchmark_bootstrap_curve)) vbm_benchmark_bootstrap_curve <- data.frame()
if (is.null(msm_benchmark_bootstrap_curve)) msm_benchmark_bootstrap_curve <- data.frame()

rgm_sharp_output <- find_rgm_sharp_T_bootstrap(
  data = analysis,
  config = config
)

rgm_sharp_bootstrap_curve <- rgm_sharp_output$results
rgm_sharp_T_star <- rgm_sharp_output$T_star

sharp_det_grid <- sort(unique(rgm_sharp_bootstrap_curve$T[is.finite(rgm_sharp_bootstrap_curve$T)]))
rgm_sharp_deterministic_curve <- run_rgm_sharp_att_T_grid(
  analysis = analysis,
  T_grid = sharp_det_grid,
  config = config
)

rgm_benchmark_output <- run_rgm_benchmark_comparison(
  analysis = analysis,
  tau_hat = tau_hat,
  full_formula = ps_formula,
  config = config,
  base_benchmark_df = benchmark_df,
  base_vbm_benchmark_bootstrap_curve = vbm_benchmark_bootstrap_curve,
  base_msm_benchmark_bootstrap_curve = msm_benchmark_bootstrap_curve
)

rgm_benchmark_df <- rgm_benchmark_output$benchmark_df

method_comparison <- data.frame(
  method = c("VBM", "MSM", "RGM-sharp", "RGM-conservative"),
  sensitivity_parameter = c("R2", "Gamma", "T", "T"),
  geometry = c("Variance of weights", "Bounded odds-ratio multiplier", "Total variation/L1 mass shift", "Total variation/L1 range bound"),
  main_assumption = c(
    "Hidden perturbation is measured through extra weight variance.",
    "Ideal control weights differ from observed weights by a bounded multiplier.",
    "The normalized ideal control-weight distribution lies within TV distance T.",
    "The same TV radius is combined with an outcome range bound."
  ),
  strength = c(
    "Interpretable covariate benchmark scale.",
    "Classic sharp marginal sensitivity model.",
    "Finite-sample transport interpretation: move probability mass from low to high outcomes.",
    "Simple closed-form conservative bound."
  ),
  weakness = c(
    "Can be unstable when observed weight variance is near zero.",
    "Gamma calibration can be conservative and ratio-driven.",
    "T calibration requires distributional benchmark interpretation.",
    "May be loose when the sample outcome range is wide."
  ),
  stringsAsFactors = FALSE
)

threshold_summary <- data.frame(
  quantity = c("ATT", "VBM_R2_star", "RGM_sharp_T_star"),
  value = c(tau_hat, vbm_R2_star, rgm_sharp_T_star),
  interpretation = c(
    "Observed IPW ATT under the NHANES observed PS model.",
    "Smallest VBM R2 on the searched grid whose bootstrap interval reaches the null.",
    "Smallest sharp RGM total-variation radius on the searched grid whose bootstrap interval reaches the null."
  ),
  stringsAsFactors = FALSE
)

write.csv(vbm_bootstrap_curve, file.path(dirs$tables, "plot_vbm_bootstrap_curve.csv"), row.names = FALSE)
write.csv(benchmark_df, file.path(dirs$tables, "plot_covariate_benchmark_vbm_msm.csv"), row.names = FALSE)
write.csv(vbm_benchmark_bootstrap_curve, file.path(dirs$tables, "plot_vbm_benchmark_bootstrap_curve.csv"), row.names = FALSE)
write.csv(msm_benchmark_bootstrap_curve, file.path(dirs$tables, "plot_msm_benchmark_bootstrap_curve.csv"), row.names = FALSE)
write.csv(rgm_sharp_bootstrap_curve, file.path(dirs$tables, "plot_rgm_sharp_bootstrap_curve.csv"), row.names = FALSE)
write.csv(rgm_sharp_deterministic_curve, file.path(dirs$tables, "plot_rgm_sharp_deterministic_curve.csv"), row.names = FALSE)
write.csv(rgm_benchmark_df, file.path(dirs$tables, "plot_covariate_benchmark_vbm_msm_rgm.csv"), row.names = FALSE)
write.csv(method_comparison, file.path(dirs$tables, "rgm_method_comparison.csv"), row.names = FALSE)
write.csv(threshold_summary, file.path(dirs$tables, "rgm_threshold_summary.csv"), row.names = FALSE)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  p_vbm <- plot_bootstrap_curve(vbm_bootstrap_curve, R2_star = vbm_R2_star) +
    ggplot2::labs(
      title = "Extension 04: VBM bootstrap curve",
      x = expression(R^2),
      y = "Bootstrap ATT interval"
    )
  ggplot2::ggsave(file.path(dirs$figures, "extension04_vbm_bootstrap_curve.png"), p_vbm, width = 7.5, height = 5.2, dpi = 300)

  p_rgm <- plot_rgm_sharp_att_T_curve(
    results = rgm_sharp_bootstrap_curve,
    T_star = rgm_sharp_T_star,
    deterministic_curve = rgm_sharp_deterministic_curve
  )
  ggplot2::ggsave(file.path(dirs$figures, "extension04_rgm_sharp_att_T_curve.png"), p_rgm, width = 7.5, height = 5.2, dpi = 300)

  p_benchmark <- plot_rgm_benchmark_comparison(
    df = rgm_benchmark_df,
    tau_hat = tau_hat,
    config = config
  )
  ggplot2::ggsave(file.path(dirs$figures, "extension04_covariate_benchmark_vbm_msm_rgm.png"), p_benchmark, width = 11, height = 6, dpi = 300)
}

v10_run_python_plot(
  project_dir = .project_dir,
  script_relative_path = "scripts/enhancedplot_rgm_vbm_msm.py",
  input_dir = dirs$tables,
  output_dir = dirs$enhanced,
  label = "Extension 04 enhanced plots"
)

cat("\n========== Extension 04 summary ==========\n")
print(threshold_summary)
cat("Outputs saved under: ", file.path("output", "extension_results", "04_rgm_vbm_msm"), "\n", sep = "")
