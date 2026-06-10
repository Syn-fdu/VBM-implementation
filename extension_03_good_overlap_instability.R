# ============================================================
# Version10 Extension 03:
# Good-overlap denominator instability
# ============================================================
# Research question:
#   Very good overlap can make observed ATT weights nearly constant. Because
#   VBM uses observed weight variance in a denominator, this can create
#   structural instability even when ordinary IPW overlap looks excellent.
#
# This Version10 script adds:
#   1. alpha x hidden-strength robustness grid,
#   2. sample-size robustness grid,
#   3. epsilon stabilizer sensitivity grid,
#   4. automatic enhanced-plot launch.
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
    "functions/ps_model.R",
    "functions/good_overlap_instability.R"
  ),
  env = .project_env
)

dirs <- v10_prepare_extension_dirs(.project_dir, "03_good_overlap_instability")
v10_clear_output_files(dirs$tables, dirs$figures, dirs$enhanced)

observed_strength_grid <- c(0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.40)
hidden_strength_grid <- c(0.03, 0.06, 0.12, 0.24, 0.36)
sample_size_grid <- c(300, 600, 1200, 1800, 3000)
epsilon_fraction_grid <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2)

cat("Running alpha x hidden-strength good-overlap grid...\n")
robustness_replicates <- do.call(
  rbind,
  lapply(hidden_strength_grid, function(h) {
    out <- run_good_overlap_instability_grid(
      observed_strength_grid = observed_strength_grid,
      replicates = 35,
      n = 1800,
      hidden_strength = h,
      seed = 90210 + round(10000 * h),
      epsilon_fraction = 1e-4,
      variance_floor = 1e-10
    )
    out$grid_type <- "alpha_by_hidden_strength"
    out
  })
)

summary_by_hidden <- do.call(
  rbind,
  lapply(split(robustness_replicates, robustness_replicates$hidden_strength), summarize_good_overlap_grid)
)
summary_by_hidden <- summary_by_hidden[order(summary_by_hidden$hidden_strength, summary_by_hidden$observed_strength), ]

cat("Running sample-size robustness grid...\n")
sample_size_replicates <- do.call(
  rbind,
  lapply(sample_size_grid, function(n_i) {
    out <- run_good_overlap_instability_grid(
      observed_strength_grid = c(0, 0.002, 0.01, 0.05, 0.20),
      replicates = 30,
      n = n_i,
      hidden_strength = 0.12,
      seed = 120000 + n_i,
      epsilon_fraction = 1e-4,
      variance_floor = 1e-10
    )
    out$sample_size <- n_i
    out$grid_type <- "sample_size"
    out
  })
)

summarize_by_two <- function(df, key1, key2) {
  pieces <- lapply(split(df, interaction(df[[key1]], df[[key2]], drop = TRUE)), function(d) {
    data.frame(
      key1_value = unique(d[[key1]])[1],
      key2_value = unique(d[[key2]])[1],
      n_replicates = nrow(d),
      flag_rate = mean(d$flag_good_overlap_instability, na.rm = TRUE),
      var_observed_control_median = stats::median(d$var_observed_control[is.finite(d$var_observed_control)], na.rm = TRUE),
      original_variance_ratio_median = stats::median(d$original_variance_ratio[is.finite(d$original_variance_ratio)], na.rm = TRUE),
      stabilized_variance_ratio_median = stats::median(d$stabilized_variance_ratio[is.finite(d$stabilized_variance_ratio)], na.rm = TRUE),
      ps_observed_iqr_median = stats::median(d$ps_observed_iqr[is.finite(d$ps_observed_iqr)], na.rm = TRUE),
      weight_cv_control_median = stats::median(d$weight_cv_control[is.finite(d$weight_cv_control)], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, pieces)
  names(out)[names(out) == "key1_value"] <- key1
  names(out)[names(out) == "key2_value"] <- key2
  out[order(out[[key1]], out[[key2]]), , drop = FALSE]
}

summary_by_sample_size <- summarize_by_two(sample_size_replicates, "sample_size", "observed_strength")

cat("Running epsilon stabilizer sensitivity grid...\n")
epsilon_replicates <- do.call(
  rbind,
  lapply(epsilon_fraction_grid, function(eps) {
    out <- run_good_overlap_instability_grid(
      observed_strength_grid = c(0, 0.001, 0.005, 0.02, 0.10),
      replicates = 30,
      n = 1800,
      hidden_strength = 0.12,
      seed = 220000 + round(-log10(eps) * 1000),
      epsilon_fraction = eps,
      variance_floor = 1e-10
    )
    out$epsilon_fraction <- eps
    out$grid_type <- "epsilon_fraction"
    out
  })
)

summary_by_epsilon <- summarize_by_two(epsilon_replicates, "epsilon_fraction", "observed_strength")

trigger_data <- make_good_overlap_instability_data(
  n = 1800,
  observed_strength = 0,
  hidden_strength = 0.12,
  seed = 90210
)

trigger_metrics <- good_overlap_metric_row(
  data = trigger_data,
  observed_strength = 0,
  hidden_strength = 0.12,
  epsilon_fraction = 1e-4,
  variance_floor = 1e-10
)

diagnostic_protocol <- data.frame(
  step = 1:5,
  diagnostic = c(
    "Report Var(w_obs | Z=0)",
    "Report control-weight CV",
    "Report original variance ratio Var(w_ideal)/Var(w_obs)",
    "Label VBM R2 unstable when the denominator is near zero or the ratio is extreme",
    "Use stabilized ratio only as supplementary numerical evidence"
  ),
  version10_rule = c(
    "var_observed_control",
    "weight_cv_control",
    "original_variance_ratio",
    "flag_good_overlap_instability",
    "stabilized_variance_ratio"
  ),
  stringsAsFactors = FALSE
)

write.csv(robustness_replicates, file.path(dirs$tables, "good_overlap_alpha_hidden_replicates.csv"), row.names = FALSE)
write.csv(summary_by_hidden, file.path(dirs$tables, "good_overlap_alpha_hidden_summary.csv"), row.names = FALSE)
write.csv(sample_size_replicates, file.path(dirs$tables, "good_overlap_sample_size_replicates.csv"), row.names = FALSE)
write.csv(summary_by_sample_size, file.path(dirs$tables, "good_overlap_sample_size_summary.csv"), row.names = FALSE)
write.csv(epsilon_replicates, file.path(dirs$tables, "good_overlap_epsilon_replicates.csv"), row.names = FALSE)
write.csv(summary_by_epsilon, file.path(dirs$tables, "good_overlap_epsilon_summary.csv"), row.names = FALSE)
write.csv(trigger_data, file.path(dirs$tables, "good_overlap_trigger_dataset.csv"), row.names = FALSE)
write.csv(trigger_metrics, file.path(dirs$tables, "good_overlap_trigger_summary.csv"), row.names = FALSE)
write.csv(diagnostic_protocol, file.path(dirs$tables, "good_overlap_diagnostic_protocol.csv"), row.names = FALSE)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  heat_df <- summary_by_hidden
  heat_df$observed_strength_f <- factor(as.character(heat_df$observed_strength), levels = as.character(observed_strength_grid))
  heat_df$hidden_strength_f <- factor(as.character(heat_df$hidden_strength), levels = as.character(hidden_strength_grid))

  p_heat <- ggplot2::ggplot(
    heat_df,
    ggplot2::aes(x = observed_strength_f, y = hidden_strength_f, fill = flag_rate)
  ) +
    ggplot2::geom_tile() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Good-overlap instability flag rate",
      subtitle = "Rows vary hidden strength; columns vary observed propensity strength alpha.",
      x = "Observed propensity strength alpha",
      y = "Hidden strength beta_U",
      fill = "Flag rate"
    )

  ggplot2::ggsave(file.path(dirs$figures, "good_overlap_alpha_hidden_flag_heatmap.png"), p_heat, width = 9.5, height = 5.8, dpi = 300)

  p_var <- ggplot2::ggplot(
    heat_df,
    ggplot2::aes(x = observed_strength, y = var_observed_control_median, group = hidden_strength_f)
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::scale_y_log10() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Observed weight variance collapses under very good overlap",
      x = "Observed propensity strength alpha",
      y = "Median Var(w_obs | Z=0), log scale"
    )

  ggplot2::ggsave(file.path(dirs$figures, "good_overlap_denominator_collapse.png"), p_var, width = 8.5, height = 5.5, dpi = 300)

  p_n <- ggplot2::ggplot(
    summary_by_sample_size,
    ggplot2::aes(x = observed_strength, y = flag_rate, group = factor(sample_size))
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Sample size does not remove structural denominator risk",
      x = "Observed propensity strength alpha",
      y = "Instability flag rate",
      group = "n"
    )

  ggplot2::ggsave(file.path(dirs$figures, "good_overlap_sample_size_robustness.png"), p_n, width = 8.5, height = 5.5, dpi = 300)
}

v10_run_python_plot(
  project_dir = .project_dir,
  script_relative_path = "scripts/enhancedplot_good_overlap_instability.py",
  input_dir = dirs$tables,
  output_dir = dirs$enhanced,
  label = "Extension 03 enhanced plots"
)

cat("\n========== Extension 03 summary ==========\n")
print(head(summary_by_hidden, 12))
cat("Outputs saved under: ", file.path("output", "version10_extension_results", "03_good_overlap_instability"), "\n", sep = "")
