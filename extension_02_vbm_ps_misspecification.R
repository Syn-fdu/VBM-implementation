# ============================================================
# Version10 Extension 02:
# VBM under observed propensity-score model misspecification
# ============================================================
# Standalone run:
#   source("extension_02_vbm_ps_misspecification.R")
#
# Research design:
#   three targeted simulation scenarios:
#     1) no confounding + correctly specified PS model;
#     2) omitted confounding strong enough to threaten the ATT conclusion;
#     3) no hidden confounding, but deliberately misspecified observed PS model.
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

# Generate one simulated dataset.
#
# Revised design for the paper claim:
#   1) no_disturbance:
#      no hidden confounder enters treatment or outcome, and the fitted
#      observed PS model is the true observed PS model.  The treatment effect
#      is deliberately strong and the outcome noise is modest, so the VBM
#      threshold should be high: VBM operates normally and the result is robust.
#
#   2) confounding_only:
#      an omitted U affects both treatment and the untreated potential outcome
#      strongly enough to move the observed ATT toward the null.  The analysis
#      PS model is correctly specified for observed covariates, but necessarily
#      omits U.  VBM should therefore report that a VBM-scale disturbance large
#      enough to overturn the conclusion is plausible/active.
#
#   3) model_misspecification_only:
#      there is no hidden confounder.  The true treatment mechanism contains
#      observed nonlinearities/interactions that are correlated with the
#      outcome, but the analysis intentionally fits a linear PS model.  This
#      lets model error produce a VBM response close to the confounding row,
#      illustrating that VBM measures perturbation strength, not its source.
generate_scenario_data <- function(
    n = 1400,
    seed = 4210,
    scenario = c("no_disturbance", "confounding_only", "model_misspecification_only"),
    hidden_treatment_strength = 1.20,
    hidden_outcome_strength = -1.05,
    misspec_treatment_strength = 1.65,
    misspec_outcome_strength = -1.35,
    true_att = 1.20,
    sigma_y = 0.55
) {
  scenario <- match.arg(scenario)
  set.seed(seed)

  gender <- factor(sample(c("Female", "Male"), n, replace = TRUE, prob = c(0.54, 0.46)))
  race <- factor(sample(c("White", "Black", "Mexican American", "Other"), n, replace = TRUE,
                        prob = c(0.55, 0.18, 0.17, 0.10)))
  age_z <- opm_std(stats::rnorm(n, mean = 45, sd = 15))
  age_z2 <- age_z^2
  income_z <- opm_std(0.30 * age_z + ifelse(gender == "Male", 0.12, -0.08) + stats::rnorm(n))
  smoking_z <- opm_std(0.35 * (gender == "Male") - 0.15 * income_z + stats::rnorm(n))

  # Always generated for diagnostics/oracle comparisons.  Its coefficients
  # are exactly zero except in the confounding_only scenario.
  U <- opm_std(
    0.25 * age_z -
      0.20 * income_z +
      ifelse(race == "Mexican American", 0.25, 0) +
      stats::rnorm(n)
  )

  eta_linear <- 0.62 * age_z +
    0.44 * income_z +
    0.28 * smoking_z +
    0.30 * (gender == "Male") -
    0.28 * (race == "Black") +
    0.24 * (race == "Mexican American")

  # Nonlinear observed PS structure.  This is present only in the model-error
  # scenario.  The correct observed formula includes these terms; the analysis
  # formula used in that row intentionally omits them.
  # H is an observed but deliberately omitted nonlinear/interaction score.
  # In the model-misspecification row, H plays the same practical role as U
  # does in the confounding row: it drives both treatment assignment and Y0.
  # The crucial difference is that H is constructed entirely from observed X,
  # so the apparent VBM disturbance is caused by a wrong PS formula, not by
  # hidden confounding.
  H <- opm_std(
    1.05 * age_z2 +
      0.95 * age_z * income_z -
      0.85 * smoking_z * (gender == "Female") +
      0.75 * income_z * (race == "Mexican American") -
      0.65 * age_z * smoking_z
  )

  eta_nonlinear <- misspec_treatment_strength * H

  eta_observed <- eta_linear
  if (identical(scenario, "model_misspecification_only")) {
    eta_observed <- eta_observed + eta_nonlinear
  }

  u_trt <- if (identical(scenario, "confounding_only")) hidden_treatment_strength else 0
  u_y <- if (identical(scenario, "confounding_only")) hidden_outcome_strength else 0

  eta_full <- eta_observed + u_trt * U
  intercept <- opm_intercept(eta_full, target_rate = 0.25)
  p_treat <- opm_inv_logit(intercept + eta_full)
  Z <- stats::rbinom(n, 1, p_treat)

  # Baseline Y0 is mostly linear.  The same active score that is omitted from
  # the analysis PS is added only in its corresponding active-source scenario:
  #   - confounding_only: hidden U affects Y0;
  #   - model_misspecification_only: observed nonlinear H affects Y0.
  # This makes the two rows look similar to VBM while preserving different
  # ground-truth causes.
  mu0_observed <- 0.42 * age_z +
    0.36 * income_z -
    0.24 * smoking_z +
    0.38 * (gender == "Female") +
    0.24 * (race == "Black") -
    0.18 * (race == "Mexican American")

  h_y <- if (identical(scenario, "model_misspecification_only")) misspec_outcome_strength else 0
  mu0 <- mu0_observed + u_y * U + h_y * H

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
    age_income = as.numeric(age_z * income_z),
    age_smoking = as.numeric(age_z * smoking_z),
    hidden_confounder = as.numeric(U),
    active_misspec_score = as.numeric(H),
    Y0 = as.numeric(Y0),
    Y1 = as.numeric(Y1),
    true_unit_ATT = as.numeric(Y1 - Y0),
    scenario_dgp = scenario,
    stringsAsFactors = FALSE
  )
}

correct_observed_ps_formula <- Z ~ gender + race + age_z + age_z2 +
  income_z * race + smoking_z * gender + age_z:income_z + age_z:smoking_z

# Deliberately simple observed PS.  This formula is correct for the no-disturbance
# and confounding-only rows because those DGPs have no observed nonlinear PS
# component.  It is intentionally wrong for model_misspecification_only.
misspecified_linear_ps_formula <- Z ~ gender + race + age_z + income_z + smoking_z

correct_oracle_full_formula <- update(correct_observed_ps_formula, . ~ . + hidden_confounder)

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

compute_misspec_target_R2 <- function(seed, config) {
  dat <- generate_scenario_data(n = 1400, seed = seed, scenario = "model_misspecification_only")
  compute_R2_between_formulas(
    full_formula = correct_observed_ps_formula,
    reduced_formula = misspecified_linear_ps_formula,
    data = dat,
    config = config
  )
}

calibrate_hidden_strength <- function(target_R2, seed, config) {
  if (!is.finite(target_R2)) return(0.90)

  candidate_strengths <- seq(0.40, 2.40, by = 0.05)
  candidate_R2 <- vapply(candidate_strengths, function(s) {
    dat <- generate_scenario_data(
      n = 1400,
      seed = seed,
      scenario = "confounding_only",
      hidden_treatment_strength = s,
      hidden_outcome_strength = -1.05
    )
    compute_R2_between_formulas(
      full_formula = correct_oracle_full_formula,
      reduced_formula = correct_observed_ps_formula,
      data = dat,
      config = config
    )
  }, numeric(1))

  best <- which.min(abs(candidate_R2 - target_R2))
  candidate_strengths[best]
}

compute_observed_minus_oracle_shift <- function(dat, analysis_formula, config) {
  ps_fit <- fit_ps_model(formula = analysis_formula, data = dat, config = config)
  ps <- clip_ps(ps_fit$fitted.values, config$ps_clip_lower, config$ps_clip_upper)
  w <- compute_weights(ps = ps, Z = dat$Z, config = config)
  tau_hat <- estimate_att(dat$Y, dat$Z, w)
  oracle_att <- estimate_from_formula(dat, correct_oracle_full_formula, config)
  tau_hat - oracle_att
}

calibrate_hidden_outcome_strength <- function(target_shift, seed, hidden_strength, config) {
  if (!is.finite(target_shift)) return(-1.05)

  candidate_y <- seq(-0.35, -2.20, by = -0.05)
  candidate_shift <- vapply(candidate_y, function(sy) {
    dat <- generate_scenario_data(
      n = 1400,
      seed = seed,
      scenario = "confounding_only",
      hidden_treatment_strength = hidden_strength,
      hidden_outcome_strength = sy
    )
    compute_observed_minus_oracle_shift(dat, correct_observed_ps_formula, config)
  }, numeric(1))

  best <- which.min(abs(candidate_shift - target_shift))
  candidate_y[best]
}

run_scenario_cell <- function(scenario_id, scenario_label, hidden_present, ps_setting, seed, hidden_strength, hidden_outcome_strength) {
  analysis <- generate_scenario_data(
    n = 1400,
    seed = seed,
    scenario = scenario_id,
    hidden_treatment_strength = hidden_strength,
    hidden_outcome_strength = hidden_outcome_strength
  )

  analysis_formula <- if (identical(ps_setting, "correct")) {
    correct_observed_ps_formula
  } else {
    misspecified_linear_ps_formula
  }

  # find_R2_bootstrap uses the global ps_formula, so set it only inside the
  # extension cell being evaluated and restore no other global state.
  ps_formula <<- analysis_formula

  ps_fit <- fit_ps_model(formula = analysis_formula, data = analysis, config = config)
  analysis$ps <- clip_ps(ps_fit$fitted.values, config$ps_clip_lower, config$ps_clip_upper)
  analysis$w <- compute_weights(ps = analysis$ps, Z = analysis$Z, config = config)

  tau_hat <- estimate_att(analysis$Y, analysis$Z, analysis$w)
  oracle_att <- estimate_from_formula(analysis, correct_oracle_full_formula, config)
  true_att <- mean(analysis$true_unit_ATT[analysis$Z == 1])

  boot <- find_R2_bootstrap(data = analysis, config = config)
  curves <- boot$results
  curves$scenario <- scenario_id
  curves$scenario_label <- scenario_label

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

  active_error_R2 <- if (!hidden_present && identical(ps_setting, "correct")) {
    0
  } else if (!hidden_present && identical(ps_setting, "misspecified")) {
    misspec_only_R2
  } else {
    hidden_only_R2
  }

  response <- if (identical(scenario_id, "no_disturbance")) {
    if (!is.finite(boot$R2_star) || boot$R2_star >= 0.25) {
      "Correctly specified/no-confounding design: VBM remains stable and the conclusion is robust"
    } else {
      "Correctly specified/no-confounding design: VBM finds only a limited sensitivity threshold"
    }
  } else if (is.finite(boot$R2_star) && is.finite(active_error_R2) && active_error_R2 >= boot$R2_star) {
    "VBM flags an active disturbance large enough to overturn the observed estimate"
  } else if (is.finite(active_error_R2) && active_error_R2 >= 0.10) {
    "VBM detects a substantial active weight-distribution disturbance"
  } else {
    "VBM detects little active disturbance on this run"
  }

  summary <- data.frame(
    scenario = scenario_id,
    scenario_label = scenario_label,
    hidden_present = hidden_present,
    ps_setting = ps_setting,
    hidden_treatment_strength = hidden_strength,
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
    active_error_R2 = active_error_R2,
    VBM_response = response,
    source_identifiable_from_VBM = "No: VBM measures weight-distribution perturbation, not its causal source",
    stringsAsFactors = FALSE
  )

  list(summary = summary, curves = curves, analysis = analysis)
}

# Three datasets requested for this extension.
misspec_seed <- config$seed + 200
hidden_seed <- config$seed + 300
no_disturbance_seed <- config$seed + 100

misspec_target_R2 <- compute_misspec_target_R2(misspec_seed, config)
hidden_strength <- calibrate_hidden_strength(misspec_target_R2, hidden_seed, config)

misspec_calibration_data <- generate_scenario_data(
  n = 1400,
  seed = misspec_seed,
  scenario = "model_misspecification_only"
)
misspec_target_shift <- compute_observed_minus_oracle_shift(
  misspec_calibration_data,
  misspecified_linear_ps_formula,
  config
)
hidden_outcome_strength <- calibrate_hidden_outcome_strength(
  target_shift = misspec_target_shift,
  seed = hidden_seed,
  hidden_strength = hidden_strength,
  config = config
)

design <- data.frame(
  scenario = c("no_disturbance", "model_misspecification_only", "confounding_only"),
  scenario_label = c("No disturbance", "Model misspecification only", "Confounding only"),
  hidden_present = c(FALSE, FALSE, TRUE),
  ps_setting = c("correct", "misspecified", "correct"),
  seed = c(no_disturbance_seed, misspec_seed, hidden_seed),
  hidden_strength = c(0, 0, hidden_strength),
  hidden_outcome_strength = c(0, 0, hidden_outcome_strength),
  stringsAsFactors = FALSE
)

runs <- lapply(seq_len(nrow(design)), function(i) {
  run_scenario_cell(
    scenario_id = design$scenario[i],
    scenario_label = design$scenario_label[i],
    hidden_present = design$hidden_present[i],
    ps_setting = design$ps_setting[i],
    seed = design$seed[i],
    hidden_strength = design$hidden_strength[i],
    hidden_outcome_strength = design$hidden_outcome_strength[i]
  )
})

summary_df <- do.call(rbind, lapply(runs, function(x) x$summary))
curves_df <- do.call(rbind, lapply(runs, function(x) x$curves))

baseline_error <- summary_df$active_error_R2[summary_df$scenario == "no_disturbance"][1]
misspec_error <- summary_df$active_error_R2[summary_df$scenario == "model_misspecification_only"][1]
confounding_error <- summary_df$active_error_R2[summary_df$scenario == "confounding_only"][1]
mimic_gap <- abs(misspec_error - confounding_error)
mimic_tolerance <- 0.05

decomposition_df <- summary_df
decomposition_df$R2_above_baseline <- decomposition_df$active_error_R2 - baseline_error
decomposition_df$error_source_label <- decomposition_df$scenario_label
decomposition_df$model_misspec_mimics_confounding <- is.finite(mimic_gap) && mimic_gap <= mimic_tolerance
decomposition_df$mimic_gap <- mimic_gap
decomposition_df$mimic_tolerance <- mimic_tolerance
decomposition_df$extension_conclusion <- ifelse(
  decomposition_df$model_misspec_mimics_confounding,
  "Model misspecification alone can generate a VBM response close to the response under confounding; VBM should not be used to identify the disturbance source.",
  "The two active-source rows are not numerically close in this random run, but VBM still reports weight-distribution perturbation rather than identifying whether it came from misspecification or confounding."
)

vbm_response_table <- decomposition_df[, c(
  "scenario",
  "scenario_label",
  "active_error_R2",
  "VBM_R2_star",
  "observed_analysis_ATT",
  "oracle_full_ATT",
  "observed_minus_oracle_ATT",
  "VBM_response",
  "source_identifiable_from_VBM",
  "extension_conclusion"
)]

write.csv(summary_df, file.path(dirs$tables, "vbm_ps_misspec_factorial_summary.csv"), row.names = FALSE)
write.csv(decomposition_df, file.path(dirs$tables, "vbm_ps_misspec_decomposition.csv"), row.names = FALSE)
write.csv(curves_df, file.path(dirs$tables, "vbm_ps_misspec_bootstrap_curves.csv"), row.names = FALSE)
write.csv(vbm_response_table, file.path(dirs$tables, "vbm_ps_misspec_vbm_response_table.csv"), row.names = FALSE)

dataset_dir <- file.path(dirs$tables, "generated_datasets")
dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)
for (i in seq_along(runs)) {
  write.csv(
    runs[[i]]$analysis,
    file.path(dataset_dir, paste0(design$scenario[i], ".csv")),
    row.names = FALSE
  )
}

if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_summary <- summary_df
  plot_summary$cell <- factor(
    plot_summary$scenario,
    levels = design$scenario,
    labels = design$scenario_label
  )

  p_r2 <- ggplot2::ggplot(plot_summary, ggplot2::aes(x = cell, y = active_error_R2)) +
    ggplot2::geom_col() +
    ggplot2::geom_point(ggplot2::aes(y = VBM_R2_star), size = 3) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1)) +
    ggplot2::labs(
      title = "VBM response under three disturbance scenarios",
      subtitle = "Bars: active source R2 on VBM scale; points: VBM R2* threshold",
      x = "Generated dataset",
      y = "R2 scale"
    )

  ggplot2::ggsave(file.path(dirs$figures, "vbm_ps_misspec_r2_decomposition.png"), p_r2, width = 9, height = 5.5, dpi = 300)

  p_att <- ggplot2::ggplot(plot_summary, ggplot2::aes(x = cell)) +
    ggplot2::geom_point(ggplot2::aes(y = observed_analysis_ATT, shape = "Observed analysis"), size = 3) +
    ggplot2::geom_point(ggplot2::aes(y = oracle_full_ATT, shape = "Oracle full"), size = 3) +
    ggplot2::geom_point(ggplot2::aes(y = true_ATT_among_treated, shape = "True ATT"), size = 3) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1)) +
    ggplot2::labs(
      title = "ATT distortion under three disturbance scenarios",
      x = "Generated dataset",
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
cat("\n========== Extension 02 VBM response table ==========\n")
print(vbm_response_table)
cat("Outputs saved under: ", file.path("output", "extension_results", "02_vbm_ps_misspecification"), "\n", sep = "")
