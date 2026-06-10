# ============================================================
# Omitted-confounder sensitivity under PS model misspecification
# Current project PS-misspecification extension
# ============================================================
#
# This file is used only by:
#
#   omitted_ps_misspec_extension/main_omitted_ps_misspec_vbm.R
#
# It does not change main.R, the NHANES reproduction workflow, the RGM
# extension, or the synthetic omitted-confounding extension.
#
# Research question:
#   When the propensity score model for the observed covariates is
#   misspecified, does VBM's benchmark logic for an omitted confounder
#   become distorted enough to change the conclusion?
#
# The omitted confounder remains the main target.  PS misspecification is
# treated as a nuisance perturbation that may interfere with the omitted-U
# benchmark.
# ============================================================

opm_inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

opm_standardize <- function(x) {
  z <- as.numeric(scale(x))
  z[!is.finite(z)] <- 0
  z
}

opm_solve_logit_intercept <- function(eta, target_rate = 0.25) {
  f <- function(a) mean(opm_inv_logit(a + eta)) - target_rate
  stats::uniroot(f, interval = c(-20, 20))$root
}

omitted_ps_misspec_design <- function() {
  data.frame(
    dataset_id = c(
      "opm_s1_mild_observed_misspec",
      "opm_s2_strong_observed_misspec",
      "opm_s3_severe_observed_misspec"
    ),
    dataset_label = c(
      "S1: mild observed-PS misspecification",
      "S2: strong observed-PS misspecification",
      "S3: severe observed-PS misspecification"
    ),
    scenario_order = 1:3,
    observed_misspecification_strength = c(0.50, 1.50, 2.00),
    hidden_u_treatment_strength = c(0.30, 0.30, 0.30),
    hidden_u_outcome_strength = c(0.25, 0.25, 0.25),
    true_att = c(0.80, 0.80, 0.80),
    sigma_y = c(1.10, 1.10, 1.10),
    seed = c(50, 50, 50),
    stringsAsFactors = FALSE
  )
}

generate_omitted_ps_misspec_profile <- function(profile, n = 650) {

  set.seed(profile$seed)

  gender_num <- stats::rbinom(n, 1, 0.48)
  gender <- factor(ifelse(gender_num == 1, "Male", "Female"))

  age_raw <- pmin(80, pmax(18, stats::rnorm(n, mean = 45, sd = 15)))
  age_z <- opm_standardize(age_raw)

  income_raw <- 45000 + 12000 * age_z + 7000 * gender_num +
    stats::rnorm(n, mean = 0, sd = 15000)
  income_raw <- pmin(120000, pmax(8000, income_raw))
  income_z <- opm_standardize(income_raw)

  race_raw <- sample(
    c("White", "Black", "Other"),
    size = n,
    replace = TRUE,
    prob = c(0.55, 0.25, 0.20)
  )
  race <- factor(race_raw, levels = c("White", "Black", "Other"))
  race_black <- as.numeric(race == "Black")
  race_other <- as.numeric(race == "Other")

  smoking_count <- stats::rpois(n, lambda = 3 + 2 * stats::rbinom(n, 1, 0.40))
  smoking_z <- opm_standardize(smoking_count)

  hidden_confounder <- opm_standardize(
    0.25 * age_z +
      0.20 * race_other -
      0.15 * race_black +
      0.20 * gender_num +
      stats::rnorm(n, mean = 0, sd = 1)
  )

  observed_linear_eta <- 0.25 * gender_num +
    0.22 * age_z +
    0.17 * income_z +
    0.12 * smoking_z +
    0.15 * race_other -
    0.12 * race_black

  observed_nonlinear_eta <- profile$observed_misspecification_strength * (
    0.45 * (age_z^2 - 1) +
      0.38 * income_z * race_other -
      0.30 * income_z * race_black +
      0.35 * smoking_z * gender_num
  )

  eta <- observed_linear_eta +
    observed_nonlinear_eta +
    profile$hidden_u_treatment_strength * hidden_confounder

  intercept <- opm_solve_logit_intercept(eta, target_rate = 0.25)
  p_treat <- opm_inv_logit(intercept + eta)
  Z <- stats::rbinom(n, size = 1, prob = p_treat)

  mu0 <- 2.50 +
    0.20 * age_z +
    0.18 * income_z +
    0.10 * smoking_z +
    0.10 * gender_num +
    0.15 * race_other -
    0.08 * race_black +
    0.18 * (age_z^2 - 1) +
    0.12 * income_z * race_other +
    profile$hidden_u_outcome_strength * hidden_confounder

  Y0 <- mu0 + stats::rnorm(n, mean = 0, sd = profile$sigma_y)
  Y1 <- Y0 + profile$true_att
  Y <- ifelse(Z == 1, Y1, Y0)

  out <- data.frame(
    Z = as.integer(Z),
    Y = as.numeric(Y),
    gender = gender,
    age_z = as.numeric(age_z),
    age_z2 = as.numeric(age_z^2),
    income_z = as.numeric(income_z),
    race = race,
    smoking_z = as.numeric(smoking_z),
    hidden_confounder = as.numeric(hidden_confounder),
    Y0 = as.numeric(Y0),
    Y1 = as.numeric(Y1),
    true_unit_ATT = as.numeric(Y1 - Y0),
    dataset_id = profile$dataset_id,
    dataset_label = profile$dataset_label,
    scenario_order = profile$scenario_order,
    observed_misspecification_strength = profile$observed_misspecification_strength,
    hidden_u_treatment_strength = profile$hidden_u_treatment_strength,
    hidden_u_outcome_strength = profile$hidden_u_outcome_strength,
    stringsAsFactors = FALSE
  )

  attr(out, "synthetic_params") <- as.list(profile)
  out
}

generate_all_omitted_ps_misspec_datasets <- function(n = 650) {
  design <- omitted_ps_misspec_design()

  out <- lapply(seq_len(nrow(design)), function(i) {
    generate_omitted_ps_misspec_profile(design[i, , drop = FALSE], n = n)
  })

  names(out) <- design$dataset_id
  attr(out, "omitted_ps_misspec_design") <- design
  out
}
