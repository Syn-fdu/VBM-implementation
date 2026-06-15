# ============================================================
# Global configuration
# NHANES VBM/MSM workflow settings with optional RGM extension
# ============================================================

config <- list(
  # Paths
  data_path = "data/nhanes.fish.rda",
  output_figures_dir = "output/figures",
  output_tables_dir = "output/tables",
  clear_output_before_run = TRUE,

  # Main-program enhanced plotting.
  # Enabled by default so source("main.R") automatically generates
  # polished enhanced versions of the base NHANES figures.
  run_main_enhanced_plots = TRUE,
  main_enhanced_output_dir = "output/figures_enhanced",

  # Reproducibility
  seed = 2024,

  # Propensity score and main ATT weights
  ps_clip_lower = 0.01,
  ps_clip_upper = 0.99,
  weight_truncation = 0.996,

  # Propensity-score covariate encoding controls
  # Numeric-like high-cardinality fields are parsed as numeric when possible;
  # categorical fields remain factors. Encoding diagnostics are printed by
  # main.R so the fitted PS model can be audited without writing extra files.
  ps_numeric_parse_min_fraction = 0.8,
  ps_high_cardinality_threshold = 15,
  ps_max_factor_levels = 10,
  ps_force_numeric_terms = c("age", "income", "income.missing", "smoking.now"),
  ps_force_factor_terms = c("gender", "race", "education", "smoking.ever"),
  ps_prefer_numeric_when_parseable = c("age", "income", "income.missing", "smoking.now"),
  ps_collapse_high_cardinality_factors = TRUE,

  # Variance-based sensitivity grid
  R2_min = 0,
  R2_max = 0.95,
  R2_coarse_step = 0.05,
  R2_fine_step = 0.01,
  R2_fine_window = 0.05,

  # Bootstrap settings for the main VBM R2* search.  One bootstrap pool is
  # generated per search stage and reused across all R2 values in that stage.
  coarse_B = 500,
  fine_B = 1000,
  percentile = c(0.025, 0.975),
  min_bootstrap_valid_fraction = 0.8,
  min_bootstrap_valid = 30,
  bootstrap_resample_by_stage = TRUE,
  diagnostics = TRUE,


  # Covariate benchmark controls
  # Benchmark propensity scores use a wider clipping range and no additional
  # benchmark-specific weight truncation so the full-vs-reduced variance and
  # Gamma calculations reflect the fitted models directly.
  benchmark_ps_clip_lower = 0.01,
  benchmark_ps_clip_upper = 0.99,
  benchmark_weight_truncation = NA_real_,
  benchmark_group = "control_weight_function",
  benchmark_min_var = 1e-10,

  # Delete observed benchmark covariates as blocks, not just single columns.
  benchmark_groups = list(
    gender = c("gender"),
    age = c("age"),
    income = c("income"),
    income_missing = c("income.missing"),
    education = c("education"),
    cig_smoked = c("smoking.now"),
    smoking_history = c("smoking.ever"),
    race = c("race")
  ),

  # Paper Figure 3-style labels and ordering for the benchmark figure.
  benchmark_plot_labels = c(
    gender = "Gender",
    age = "Age",
    income = "Income",
    income_missing = "Income (Missing)",
    education = "Education",
    cig_smoked = "Cig. Smoked",
    smoking_history = "Smoking history",
    race = "Race"
  ),

  # VBM benchmark bootstrap.  Each covariate benchmark R2_j is computed once
  # on the original sample and held fixed inside a shared bootstrap pool.
  vbm_benchmark_bootstrap_enabled = TRUE,
  vbm_benchmark_bootstrap_B = 1000,
  vbm_benchmark_bootstrap_seed_offset = 8100,

  # Paper Figure 3 benchmark extensions. These flags add rows to the
  # benchmark table without changing the base VBM/MSM parameter settings.
  figure3_include_vbm_corr = TRUE,
  figure3_include_qbal = TRUE,
  vbm_corr_group = "control_weight_function",
  msm_qbal_ratio_quantile = 0.95,

  # MSM bound engine and benchmark-bootstrap defaults
  Gamma_min = 1,
  Gamma_max = 3,
  Gamma_step = 0.05,
  msm_min_control_n = 30,
  msm_bootstrap_B = 500,

  # MSM benchmark Gamma is the raw maximum symmetric full-vs-reduced ATT
  # weight-function ratio among controls.
  msm_gamma_group = "control_weight_function",

  # MSM benchmark bootstrap.  Each covariate benchmark Gamma_j is computed
  # once on the original sample and held fixed inside a shared bootstrap pool.
  msm_benchmark_bootstrap_enabled = TRUE,
  msm_benchmark_bootstrap_B = 1000,
  msm_benchmark_bootstrap_seed_offset = 9100,

  # Fallback single-term benchmark list used only when benchmark_groups is unavailable.
  benchmark_vars = c(
    "gender",
    "age",
    "income",
    "income.missing",
    "education",
    "smoking.ever",
    "smoking.now",
    "race"
  ),

  # current project extension runner controls.
  # All four refactored extensions are enabled by default, so source("main.R")
  # runs the base NHANES workflow and then the four extension workflows.
  # Any individual extension can be disabled by editing the corresponding flag below.
  run_extensions = TRUE,
  run_extension_01_hidden_strength_vbm_msm = TRUE,
  run_extension_02_vbm_ps_misspecification = TRUE,
  run_extension_03_good_overlap_instability = TRUE,
  run_extension_04_rgm_vbm_msm = TRUE,


  # RGM sharp ATT-T search grid and bootstrap settings.
  # The sharp model is the only RGM model that searches T*.
  rgm_T_min = 0,
  rgm_T_max = 1,
  rgm_T_coarse_step = 0.05,
  rgm_T_fine_step = 0.005,
  rgm_T_fine_window = 0.05,
  rgm_sharp_coarse_B = 500,
  rgm_sharp_fine_B = 1000,
  rgm_sharp_search_seed_offset = 12100,

  # RGM conservative range-bound settings.
  # This version does not search T*; it only contributes fixed-T benchmark rows.
  rgm_conservative_outcome_bound_mode = "sample_control",
  rgm_winsor_probs = c(0.01, 0.99),
  rgm_y_lower = NA_real_,
  rgm_y_upper = NA_real_,

  # RGM benchmark settings.  T_j is computed as the total-variation distance
  # between normalized full and reduced ATT weight-function distributions.
  rgm_benchmark_bootstrap_enabled = TRUE,
  rgm_benchmark_bootstrap_B = 1000,
  rgm_benchmark_group = "control_weight_function",
  rgm_sharp_benchmark_seed_offset = 13100,
  rgm_conservative_benchmark_seed_offset = 11100

)
