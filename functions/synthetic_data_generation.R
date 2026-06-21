# ============================================================
# Synthetic NHANES-fish-like data generation for VBM extension
# Current project: ordered omitted-confounding stress test
# ============================================================
#
# This file is used only by extension_01_hidden_strength_vbm_msm.R.
# The generated datasets form an ordered omitted-confounding stress test,
# ranging from no hidden confounding to hidden confounding strong enough
# to plausibly overturn the observed-only conclusion.
# ============================================================

inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

solve_logit_intercept <- function(eta, target_rate = 0.22) {
  f <- function(a) mean(inv_logit(a + eta)) - target_rate
  stats::uniroot(f, interval = c(-20, 20))$root
}

standardize_num <- function(x) {
  z <- as.numeric(scale(x))
  z[!is.finite(z)] <- 0
  z
}

make_base_covariates <- function(n, seed) {

  set.seed(seed)

  gender <- factor(
    sample(c("Female", "Male"), size = n, replace = TRUE, prob = c(0.53, 0.47))
  )

  age <- pmin(80, pmax(18, round(stats::rnorm(n, mean = 45, sd = 15))))

  income_missing <- stats::rbinom(n, 1, 0.11)

  income_z_raw <- 0.30 * standardize_num(age) +
    ifelse(gender == "Male", 0.08, -0.05) +
    stats::rnorm(n, 0, 1)

  income <- round(55000 + 18000 * income_z_raw)
  income <- pmin(130000, pmax(8000, income))
  income[income_missing == 1] <- round(stats::median(income, na.rm = TRUE))

  race <- factor(
    sample(
      c("White", "Black", "Mexican American", "Other"),
      size = n,
      replace = TRUE,
      prob = c(0.55, 0.18, 0.17, 0.10)
    )
  )

  education <- factor(
    sample(
      c("<HS", "HS", "Some college", "College+"),
      size = n,
      replace = TRUE,
      prob = c(0.15, 0.25, 0.32, 0.28)
    ),
    levels = c("<HS", "HS", "Some college", "College+")
  )

  smoking_ever <- factor(
    sample(c("Never", "Ever"), size = n, replace = TRUE, prob = c(0.56, 0.44))
  )

  smoking_now <- ifelse(
    smoking_ever == "Ever",
    stats::rpois(n, lambda = 5.5),
    stats::rbinom(n, 1, 0.03)
  )

  location_region <- factor(
    sample(
      c("North", "South", "Coastal", "Inland"),
      size = n,
      replace = TRUE,
      prob = c(0.26, 0.25, 0.25, 0.24)
    )
  )

  data.frame(
    gender = gender,
    age = age,
    income = income,
    income.missing = income_missing,
    race = race,
    education = education,
    smoking.ever = smoking_ever,
    smoking.now = smoking_now,
    location_region = location_region,
    stringsAsFactors = FALSE
  )
}

observed_treatment_eta <- function(d, scale = 1) {

  age_z <- standardize_num(d$age)
  income_z <- standardize_num(d$income)
  smoke_z <- standardize_num(d$smoking.now)

  eta <- 0.22 * (d$gender == "Male") +
    0.25 * age_z +
    0.18 * income_z +
    0.10 * smoke_z +
    0.16 * (d$smoking.ever == "Ever") +
    0.16 * (d$education == "College+") -
    0.09 * (d$education == "<HS") +
    0.20 * (d$race == "Mexican American") -
    0.10 * (d$race == "Black") +
    0.05 * d$income.missing

  as.numeric(scale) * eta
}

observed_outcome_mean <- function(d, scale = 1) {

  age_z <- standardize_num(d$age)
  income_z <- standardize_num(d$income)
  smoke_z <- standardize_num(d$smoking.now)

  mu <- 3.00 +
    0.17 * age_z +
    0.12 * income_z +
    0.07 * smoke_z +
    0.08 * (d$gender == "Male") +
    0.12 * (d$smoking.ever == "Ever") +
    0.10 * (d$race == "Mexican American") -
    0.08 * (d$race == "Black") +
    0.08 * (d$education == "College+") -
    0.05 * (d$education == "<HS")

  3.00 + as.numeric(scale) * (mu - 3.00)
}

location_treatment_effect <- function(location_region, strength = 0) {
  vals <- c(
    North = -0.15,
    South = 0.25,
    Coastal = 0.55,
    Inland = -0.35
  )
  strength * unname(vals[as.character(location_region)])
}

location_outcome_effect <- function(location_region, strength = 0) {
  vals <- c(
    North = -0.12,
    South = 0.18,
    Coastal = 0.45,
    Inland = -0.25
  )
  strength * unname(vals[as.character(location_region)])
}

# Compact, ordered design grid. strength_index is a visual axis for
# comparing VBM and MSM detection across datasets; it is not used in the DGP.

synthetic_design_grid <- function() {
  # Project hidden-strength grid:
  # every non-baseline scenario deliberately omits the same hidden variable.
  # The only systematic change across rows is hidden treatment/outcome strength.
  data.frame(
    dataset_id = c(
      "hidden_h0_none",
      "hidden_h1_tiny",
      "hidden_h2_weak",
      "hidden_h3_moderate",
      "hidden_h4_strong",
      "hidden_h5_severe"
    ),
    dataset_label = c(
      "H0: no hidden confounding",
      "H1: tiny hidden variable",
      "H2: weak hidden variable",
      "H3: moderate hidden variable",
      "H4: strong hidden variable",
      "H5: severe hidden variable"
    ),
    scenario_class = c(
      "none",
      "tiny_hidden",
      "weak_hidden",
      "moderate_hidden",
      "strong_hidden",
      "severe_hidden"
    ),
    scenario_order = 0:5,
    strength_index = c(0.00, 0.03, 0.06, 0.12, 0.24, 0.36),
    seed = c(3101, 3201, 3301, 3401, 3501, 3601),
    omitted_variable = rep("hidden_confounder", 6),
    true_att = c(0.90, 0.90, 0.90, 0.90, 0.90, 0.90),
    sigma_y = c(0.85, 0.88, 0.92, 1.00, 1.12, 1.28),
    observed_trt_scale = c(1.00, 0.98, 0.95, 0.90, 0.82, 0.75),
    observed_y_scale = c(1.00, 0.98, 0.95, 0.90, 0.82, 0.75),
    location_trt_strength = rep(0.00, 6),
    location_y_strength = rep(0.00, 6),
    hidden_trt_strength = c(0.00, 0.25, 0.55, 1.10, 1.80, 2.60),
    hidden_y_strength = c(0.00, 0.05, 0.10, 0.20, 0.32, 0.48),
    hidden_sd = rep(1.00, 6),
    stringsAsFactors = FALSE
  )
}

generate_synthetic_analysis <- function(
    dataset_id,
    dataset_label,
    n = 540,
    seed = 2027,
    true_att = 1.15,
    sigma_y = 0.75,
    target_treated_rate = 0.22,
    observed_trt_scale = 1.00,
    observed_y_scale = 1.00,
    location_trt_strength = 0.00,
    location_y_strength = 0.00,
    hidden_trt_strength = 0.00,
    hidden_y_strength = 0.00,
    hidden_sd = 1.00,
    omitted_variable = "hidden_confounder",
    scenario_order = NA_integer_,
    scenario_class = NA_character_,
    strength_index = NA_real_
) {

  d <- make_base_covariates(n = n, seed = seed)

  set.seed(seed + 101)

  # Hidden confounder is always generated.  Its treatment/outcome coefficients
  # are zero in S0/S1 and increase monotonically from S2 to S6.
  hidden_location_mean <- c(
    North = -0.25,
    South = 0.10,
    Coastal = 0.30,
    Inland = -0.15
  )

  hidden_confounder <- as.numeric(
    hidden_location_mean[as.character(d$location_region)] +
      stats::rnorm(n, mean = 0, sd = hidden_sd)
  )

  hidden_confounder <- standardize_num(hidden_confounder)

  eta <- observed_treatment_eta(d, scale = observed_trt_scale) +
    location_treatment_effect(d$location_region, strength = location_trt_strength) +
    hidden_trt_strength * hidden_confounder

  intercept <- solve_logit_intercept(eta, target_rate = target_treated_rate)
  p_treat <- inv_logit(intercept + eta)

  Z <- stats::rbinom(n, size = 1, prob = p_treat)

  mu0 <- observed_outcome_mean(d, scale = observed_y_scale) +
    location_outcome_effect(d$location_region, strength = location_y_strength) +
    hidden_y_strength * hidden_confounder

  eps <- stats::rnorm(n, mean = 0, sd = sigma_y)

  Y0 <- mu0 + eps
  Y1 <- Y0 + true_att
  Y <- ifelse(Z == 1, Y1, Y0)

  out <- data.frame(
    Z = as.integer(Z),
    Y = as.numeric(Y),
    gender = factor(d$gender),
    age = as.numeric(d$age),
    income = as.numeric(d$income),
    income.missing = as.numeric(d$income.missing),
    race = factor(d$race),
    education = factor(d$education),
    smoking.ever = factor(d$smoking.ever),
    smoking.now = as.numeric(d$smoking.now),
    location_region = factor(d$location_region),
    hidden_confounder = as.numeric(hidden_confounder),
    Y0 = as.numeric(Y0),
    Y1 = as.numeric(Y1),
    true_unit_ATT = as.numeric(Y1 - Y0),
    dataset_id = dataset_id,
    dataset_label = dataset_label,
    scenario_order = scenario_order,
    scenario_class = scenario_class,
    omitted_confounding_strength = strength_index,
    stringsAsFactors = FALSE
  )

  omitted_terms <- omitted_variable

  attr(out, "true_att") <- true_att
  attr(out, "target_treated_rate") <- target_treated_rate
  attr(out, "p_treat") <- p_treat
  attr(out, "intercept") <- intercept
  attr(out, "omitted_variable") <- omitted_variable
  attr(out, "omitted_terms") <- omitted_terms
  attr(out, "scenario_order") <- scenario_order
  attr(out, "scenario_class") <- scenario_class
  attr(out, "omitted_confounding_strength") <- strength_index
  attr(out, "synthetic_params") <- list(
    dataset_id = dataset_id,
    dataset_label = dataset_label,
    seed = seed,
    true_att = true_att,
    sigma_y = sigma_y,
    target_treated_rate = target_treated_rate,
    observed_trt_scale = observed_trt_scale,
    observed_y_scale = observed_y_scale,
    location_trt_strength = location_trt_strength,
    location_y_strength = location_y_strength,
    hidden_trt_strength = hidden_trt_strength,
    hidden_y_strength = hidden_y_strength,
    hidden_sd = hidden_sd,
    omitted_variable = omitted_variable,
    omitted_terms = omitted_terms,
    scenario_order = scenario_order,
    scenario_class = scenario_class,
    omitted_confounding_strength = strength_index
  )

  out
}

generate_synthetic_profile <- function(profile, n = 540) {
  generate_synthetic_analysis(
    dataset_id = profile$dataset_id,
    dataset_label = profile$dataset_label,
    n = n,
    seed = profile$seed,
    true_att = profile$true_att,
    sigma_y = profile$sigma_y,
    target_treated_rate = 0.22,
    observed_trt_scale = profile$observed_trt_scale,
    observed_y_scale = profile$observed_y_scale,
    location_trt_strength = profile$location_trt_strength,
    location_y_strength = profile$location_y_strength,
    hidden_trt_strength = profile$hidden_trt_strength,
    hidden_y_strength = profile$hidden_y_strength,
    hidden_sd = profile$hidden_sd,
    omitted_variable = profile$omitted_variable,
    scenario_order = profile$scenario_order,
    scenario_class = profile$scenario_class,
    strength_index = profile$strength_index
  )
}

generate_all_synthetic_datasets <- function(n = 540) {
  design <- synthetic_design_grid()

  out <- lapply(seq_len(nrow(design)), function(i) {
    generate_synthetic_profile(design[i, , drop = FALSE], n = n)
  })

  names(out) <- design$dataset_id
  attr(out, "synthetic_design_grid") <- design
  out
}
