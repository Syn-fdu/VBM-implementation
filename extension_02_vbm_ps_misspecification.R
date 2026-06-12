# ============================================================
# Version10 Extension 02:
# VBM under observed propensity-score model misspecification
# ============================================================
# Standalone run:
#   source("extension_02_vbm_ps_misspecification.R")
#
# Research design:
#   2 x 2 factorial simulation:
#     hidden confounder: absent / present
#     observed PS model: correct / misspecified
#
# Main outputs:
#   - vbm_ps_misspec_factorial_summary.csv
#   - vbm_ps_misspec_decomposition.csv
#   - vbm_ps_misspec_bootstrap_curves.csv
#   - standard and enhanced plots
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
    "functions/vbm_bounds.R",
    "functions/bootstrap_analysis.R",
    "functions/msm_analysis.R",
    "functions/benchmark_analysis.R",
    "functions/plotting.R"
  ),
  env = .project_env
)

dirs <- v10_prepare_extension_dirs(.project_dir, "02_vbm_ps_misspecification")
v10_clear_output_files(dirs$tables, dirs$figures, dirs$enhanced)

config$seed <- 4210
config$output_figures_dir <- file.path("output", "extension_results", "02_vbm_ps_misspecification", "figures")
config$output_tables_dir <- file.path("output", "extension_results", "02_vbm_ps_misspecification", "tables")

# Keep the extension fast enough to run automatically from main.R.
config$R2_min <- 0
config$R2_max <- 0.90
config$R2_coarse_step <- 0.10
config$R2_fine_step <- 0.02
config$R2_fine_window <- 0.08
config$coarse_B <- 70
config$fine_B <- 100
config$min_bootstrap_valid_fraction <- 0.60
config$min_bootstrap_valid <- 25
config$diagnostics <- TRUE
config$weight_truncation <- 0.996
config$benchmark_weight_truncation <- NA_real_

opm_inv_logit <- function(x) 1 / (1 + exp(-pmax(pmin(x, 35), -35)))

opm_std <- function(x) {
  z <- as.numeric(scale(x))
  z[!is.finite(z)] <- 0
  z
}

opm_intercept <- function(eta, target_rate = 0.25) {
  f <- function(a) mean(opm_inv_logit(a + eta)) - target_rate
  stats::uniroot(f, c(-20, 20))$root
}

generate_factorial_data <- function(
    n = 900,
    seed = 4210,
    hidden_present = FALSE,
    hidden_treatment_strength = 0.90,
    hidden_outcome_strength = 0.40,
    true_att = 0.80,
    sigma_y = 1.10
) {
  set.seed(seed)

  gender <- factor(sample(c("Female", "Male"), n, replace = TRUE, prob = c(0.54, 0.46)))
  race <- factor(sample(c("White", "Black", "Mexican American", "Other"), n, replace = TRUE,
                        prob = c(0.55, 0.18, 0.17, 0.10)))
  age_z <- opm_std(stats::rnorm(n, mean = 45, sd = 15))
  age_z2 <- age_z^2
  income_z <- opm_std(0.30 * age_z + ifelse(gender == "Male", 0.12, -0.08) + stats::rnorm(n))
  smoking_z <- opm_std(0.35 * (gender == "Male") - 0.15 * income_z + stats::rnorm(n))

  U <- opm_std(
    0.25 * age_z -
      0.20 * income_z +
      ifelse(race == "Mexican American", 0.25, 0) +
      stats::rnorm(n)
  )

  u_trt <- if (hidden_present) hidden_treatment_strength else 0
  u_y <- if (hidden_present) hidden_outcome_strength else 0

  eta_observed <- 0.70 * age_z -
    0.40 * age_z2 +
    0.45 * income_z +
    0.35 * smoking_z +
    0.35 * (gender == "Male") -
    0.30 * (race == "Black") +
    0.25 * (race == "Mexican American") +
    0.35 * income_z * (race == "Mexican American") -
    0.30 * smoking_z * (gender == "Female")

  eta_full <- eta_observed + u_trt * U
  intercept <- opm_intercept(eta_full, target_rate = 0.25)
  p_treat <- opm_inv_logit(intercept + eta_full)
  Z <- stats::rbinom(n, 1, p_treat)

  mu0 <- 0.45 * age_z + 0.35 * income_z - 0.25 * smoking_z +
    0.40 * (gender == "Female") +
    0.25 * (race == "Black") -
    0.20 * (race == "Mexican American") +
    0.30 * age_z2 +
    u_y * U

  Y0 <- mu0 + stats::rnorm(n, sd = sigma_y)
  Y1 <- Y0 + true_att
  Y <- ifelse(Z == 1, Y1, Y0)

  data.frame(
    Y = as.numeric(Y),
    Z = as.integer(Z),
    gender = gender,
    race = race,
    age_z = as.numeric(age_z),
    age_z2 = as.numeric(age_z2),
    income_z = as.numeric(income_z),
    smoking_z = as.numeric(smoking_z),
    hidden_confounder = as.numeric(U),
    Y0 = as.numeric(Y0),
    Y1 = as.numeric(Y1),
    true_unit_ATT = as.numeric(Y1 - Y0),
    stringsAsFactors = FALSE
  )
}

correct_observed_ps_formula <- Z ~ gender + race + age_z + age_z2 +
  income_z * race + smoking_z * gender

misspecified_linear_ps_formula <- Z ~ gender + race + age_z + income_z + smoking_z

correct_oracle_full_formula <- update(correct_observed_ps_formula, . ~ . + hidden_confounder)
misspecified_plus_U_formula <- update(misspecified_linear_ps_formula, . ~ . + hidden_confounder)

estimate_from_formula <- function(analysis, formula, config) {
  fit <- fit_ps_model(formula = formula, data = analysis, config = config)
  ps <- clip_ps(fit$fitted.values, config$ps_clip_lower, config$ps_clip_upper)
  w <- compute_weights(ps = ps, Z = analysis$Z, config = config)
  estimate_att(Y = analysis$Y, Z = analysis$Z, w = w)
}

compute_R2_between_formulas <- function(full_formula, reduced_formula, data, config) {
  bcfg <- make_benchmark_config(config)

  full_fit <- fit_ps_model(full_formula, data = data, config = bcfg)
  reduced_fit <- fit_ps_model(reduced_formula, data = data, config = bcfg)

  ps_full <- clip_ps(full_fit$fitted.values, bcfg$ps_clip_lower, bcfg$ps_clip_upper)
  ps_reduced <- clip_ps(reduced_fit$fitted.values, bcfg$ps_clip_lower, bcfg$ps_clip_upper)

  w_full <- compute_att_weight_function(ps = ps_full, Z = data$Z)
  w_reduced <- compute_att_weight_function(ps = ps_reduced, Z = data$Z)

  info <- compute_benchmark_R2_value(
    w_full = w_full,
    w_reduced = w_reduced,
    data = data,
    group = "control_weight_function",
    config = bcfg
  )

  info$R2_used
}

run_factorial_cell <- function(hidden_present, ps_setting, seed) {
  analysis <- generate_factorial_data(
    n = 900,
    seed = seed,
    hidden_present = hidden_present
  )

  analysis_formula <- if (identical(ps_setting, "correct")) {
    correct_observed_ps_formula
  } else {
    misspecified_linear_ps_formula
  }

  ps_formula <<- analysis_formula

  ps_fit <- fit_ps_model(formula = analysis_formula, data = analysis, config = config)
  analysis$ps <- clip_ps(ps_fit$fitted.values, config$ps_clip_lower, config$ps_clip_upper)
  analysis$w <- compute_weights(ps = analysis$ps, Z = analysis$Z, config = config)

  tau_hat <- estimate_att(analysis$Y, analysis$Z, analysis$w)
  oracle_formula <- correct_oracle_full_formula
  oracle_att <- estimate_from_formula(analysis, oracle_formula, config)
  true_att <- mean(analysis$true_unit_ATT[analysis$Z == 1])

  boot <- find_R2_bootstrap(data = analysis, config = config)
  curves <- boot$results
  curves$scenario <- paste0(ifelse(hidden_present, "hidden", "no_hidden"), "_", ps_setting)

  hidden_only_R2 <- compute_R2_between_formulas(
    full_formula = correct_oracle_full_formula,
    reduced_formula = correct_observed_ps_formula,
    data = analysis,
    config = config
  )

  misspec_only_R2 <- compute_R2_between_formulas(
    full_formula = correct_observed_ps_formula,
    reduced_formula = misspecified_linear_ps_formula,
    data = analysis,
    config = config
  )

  combined_error_R2 <- compute_R2_between_formulas(
    full_formula = correct_oracle_full_formula,
    reduced_formula = analysis_formula,
    data = analysis,
    config = config
  )

  summary <- data.frame(
    scenario = paste0(ifelse(hidden_present, "hidden", "no_hidden"), "_", ps_setting),
    hidden_present = hidden_present,
    ps_setting = ps_setting,
    n_total = nrow(analysis),
    n_treated = sum(analysis$Z == 1),
    n_control = sum(analysis$Z == 0),
    true_ATT_among_treated = true_att,
    observed_analysis_ATT = tau_hat,
    oracle_full_ATT = oracle_att,
    observed_minus_oracle_ATT = tau_hat - oracle_att,
    VBM_R2_star = boot$R2_star,
    hidden_only_R2 = hidden_only_R2,
    misspec_only_R2 = misspec_only_R2,
    combined_error_R2 = combined_error_R2,
    interaction_gap = combined_error_R2 - hidden_only_R2 - misspec_only_R2,
    misspec_can_mimic_hidden = (!hidden_present && ps_setting == "misspecified" && is.finite(combined_error_R2) && combined_error_R2 > 0.02),
    benchmark_can_overturn = is.finite(combined_error_R2) && is.finite(boot$R2_star) && combined_error_R2 >= boot$R2_star,
    stringsAsFactors = FALSE
  )

  list(summary = summary, curves = curves, analysis = analysis)
}

design <- expand.grid(
  hidden_present = c(FALSE, TRUE),
  ps_setting = c("correct", "misspecified"),
  stringsAsFactors = FALSE
)
design$seed <- config$seed + seq_len(nrow(design)) * 100

runs <- lapply(seq_len(nrow(design)), function(i) {
  run_factorial_cell(
    hidden_present = design$hidden_present[i],
    ps_setting = design$ps_setting[i],
    seed = design$seed[i]
  )
})

summary_df <- do.call(rbind, lapply(runs, function(x) x$summary))
curves_df <- do.call(rbind, lapply(runs, function(x) x$curves))

baseline_combined <- summary_df$combined_error_R2[
  summary_df$hidden_present == FALSE & summary_df$ps_setting == "correct"
][1]

decomposition_df <- summary_df
decomposition_df$R2_above_baseline <- decomposition_df$combined_error_R2 - baseline_combined
decomposition_df$error_source_label <- ifelse(
  !decomposition_df$hidden_present & decomposition_df$ps_setting == "correct",
  "A1: neither hidden U nor misspecification",
  ifelse(
    !decomposition_df$hidden_present & decomposition_df$ps_setting == "misspecified",
    "A2: misspecification only",
    ifelse(
      decomposition_df$hidden_present & decomposition_df$ps_setting == "correct",
      "B1: hidden U only",
      "B2: hidden U + misspecification"
    )
  )
)

write.csv(summary_df, file.path(dirs$tables, "vbm_ps_misspec_factorial_summary.csv"), row.names = FALSE)
write.csv(decomposition_df, file.path(dirs$tables, "vbm_ps_misspec_decomposition.csv"), row.names = FALSE)
write.csv(curves_df, file.path(dirs$tables, "vbm_ps_misspec_bootstrap_curves.csv"), row.names = FALSE)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_summary <- summary_df
  plot_summary$cell <- factor(
    plot_summary$scenario,
    levels = c("no_hidden_correct", "no_hidden_misspecified", "hidden_correct", "hidden_misspecified")
  )

  p_r2 <- ggplot2::ggplot(plot_summary, ggplot2::aes(x = cell, y = combined_error_R2)) +
    ggplot2::geom_col() +
    ggplot2::geom_point(ggplot2::aes(y = VBM_R2_star), size = 3) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)) +
    ggplot2::labs(
      title = "VBM decomposition under PS misspecification",
      subtitle = "Bars: combined hidden/misspecification R2; points: VBM R2* threshold",
      x = "2 x 2 simulation cell",
      y = "R2 scale"
    )

  ggplot2::ggsave(file.path(dirs$figures, "vbm_ps_misspec_r2_decomposition.png"), p_r2, width = 9, height = 5.5, dpi = 300)

  p_att <- ggplot2::ggplot(plot_summary, ggplot2::aes(x = cell)) +
    ggplot2::geom_point(ggplot2::aes(y = observed_analysis_ATT, shape = "Observed analysis"), size = 3) +
    ggplot2::geom_point(ggplot2::aes(y = oracle_full_ATT, shape = "Oracle full"), size = 3) +
    ggplot2::geom_point(ggplot2::aes(y = true_ATT_among_treated, shape = "True ATT"), size = 3) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)) +
    ggplot2::labs(
      title = "ATT distortion from hidden U and PS misspecification",
      x = "2 x 2 simulation cell",
      y = "ATT",
      shape = "Quantity"
    )

  ggplot2::ggsave(file.path(dirs$figures, "vbm_ps_misspec_att_distortion.png"), p_att, width = 9, height = 5.5, dpi = 300)
}

v10_run_python_plot(
  project_dir = .project_dir,
  script_relative_path = "scripts/enhancedplot_vbm_ps_misspecification.py",
  input_dir = dirs$tables,
  output_dir = dirs$enhanced,
  label = "Extension 02 enhanced plots"
)

cat("\n========== Extension 02 summary ==========\n")
print(summary_df)
cat("Outputs saved under: ", file.path("output", "extension_results", "02_vbm_ps_misspecification"), "\n", sep = "")
