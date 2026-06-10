# ============================================================
# Propensity score model and ATT weighting
# ============================================================

ps_formula <- Z ~ gender + age + income +
  income.missing + race +
  education + smoking.ever +
  smoking.now

fit_ps_model <- function(formula, data, config = NULL) {

  X <- model.matrix(formula, data)

  maxit <- 50
  if (!is.null(config$ps_maxit)) {
    maxit <- config$ps_maxit
  }

  suppressWarnings(
    glm.fit(
      x = X,
      y = data$Z,
      family = binomial(),
      control = glm.control(maxit = maxit)
    )
  )
}

clip_ps <- function(ps, lower = 0.01, upper = 0.99) {

  pmin(
    pmax(ps, lower),
    upper
  )
}

compute_att_weight_function <- function(ps, Z) {

  # Population ATT weight function from the paper:
  #   w(X) = Pr(Z=0)/Pr(Z=1) * e(X)/(1-e(X)).
  # Unlike compute_weights(), this returns w(X) for every unit and
  # does not set treated observations to 1. It is used only for
  # covariate benchmarking and MSM Gamma diagnostics.
  p1 <- mean(Z)

  (ps / (1 - ps)) * ((1 - p1) / p1)
}

compute_weights <- function(ps, Z, config = NULL) {

  p1 <- mean(Z)

  w <- ifelse(
    Z == 0,
    (ps / (1 - ps)) * ((1 - p1) / p1),
    1
  )

  if (!is.null(config$weight_truncation)) {
    q <- config$weight_truncation

    if (!is.na(q) && q < 1) {
      upper <- quantile(
        w[Z == 0],
        probs = q,
        na.rm = TRUE,
        names = FALSE
      )

      w[Z == 0] <- pmin(
        w[Z == 0],
        upper
      )
    }
  }

  w
}

estimate_att <- function(Y, Z, w) {

  treated_mean <- mean(Y[Z == 1])

  control_mean <- weighted.mean(
    x = Y[Z == 0],
    w = w[Z == 0]
  )

  treated_mean - control_mean
}
