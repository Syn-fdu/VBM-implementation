# ============================================================
# Data preparation
# Low-dimensional PS encoding + explicit encoding diagnostics
# ============================================================

library(tidyverse)

.config_value <- function(config, name, default) {
  if (is.null(config) || is.null(config[[name]])) {
    return(default)
  }
  config[[name]]
}

missing_like_values <- c(
  "", ".", "NA", "N/A", "na", "n/a", "NaN",
  "Missing", "missing", "Unknown", "unknown",
  "Refused", "refused", "Don't know", "Don't Know",
  "dont know", "Don't know/Not sure", "Not ascertained"
)

clean_character_values <- function(x) {
  y <- as.character(x)
  y <- trimws(y)
  y[y %in% missing_like_values] <- NA_character_
  y
}

parse_numeric_like <- function(x) {

  y <- clean_character_values(x)
  y_no_comma <- gsub(",", "", y)
  y_no_currency <- gsub("[$]", "", y_no_comma)

  direct <- suppressWarnings(as.numeric(y_no_currency))
  out <- direct

  # If labels are ranges such as "20000-24999" or "$20,000 to $24,999",
  # use the midpoint of the first two numbers. If only one number appears,
  # use that number. Direct numeric strings have already been handled above.
  needs_regex <- is.na(out) & !is.na(y_no_currency)
  if (any(needs_regex)) {
    pieces <- gregexpr("\\d+\\.?\\d*", y_no_currency[needs_regex])
    nums <- regmatches(y_no_currency[needs_regex], pieces)
    idx <- which(needs_regex)

    for (k in seq_along(nums)) {
      vals <- suppressWarnings(as.numeric(nums[[k]]))
      vals <- vals[is.finite(vals)]

      if (length(vals) >= 2) {
        out[idx[k]] <- mean(vals[1:2])
      } else if (length(vals) == 1) {
        out[idx[k]] <- vals[1]
      }
    }
  }

  out
}

collapse_factor_levels <- function(x, max_levels = 10) {

  y <- clean_character_values(x)

  if (all(is.na(y))) {
    return(factor(y))
  }

  tab <- sort(table(y, useNA = "no"), decreasing = TRUE)
  keep <- names(tab)[seq_len(min(length(tab), max_levels))]

  y_collapsed <- y
  y_collapsed[!is.na(y_collapsed) & !(y_collapsed %in% keep)] <- "__other__"

  factor(y_collapsed)
}

encode_ps_variable <- function(
    x,
    variable,
    config,
    default_role = "factor"
) {

  clean <- clean_character_values(x)
  raw_class <- paste(class(x), collapse = "/")
  raw_n_unique <- length(unique(clean[!is.na(clean)]))
  n_missing_raw <- sum(is.na(clean))

  force_numeric <- variable %in% .config_value(config, "ps_force_numeric_terms", character())
  force_factor <- variable %in% .config_value(config, "ps_force_factor_terms", character())
  prefer_numeric <- variable %in% .config_value(config, "ps_prefer_numeric_when_parseable", character())

  parse_min_fraction <- .config_value(config, "ps_numeric_parse_min_fraction", 0.8)
  high_cardinality_threshold <- .config_value(config, "ps_high_cardinality_threshold", 15)
  max_factor_levels <- .config_value(config, "ps_max_factor_levels", 10)
  collapse_high_cardinality <- .config_value(config, "ps_collapse_high_cardinality_factors", TRUE)

  parsed <- parse_numeric_like(x)
  nonmissing_clean <- !is.na(clean)
  parse_fraction <- ifelse(
    sum(nonmissing_clean) > 0,
    mean(!is.na(parsed[nonmissing_clean])),
    0
  )

  use_numeric <- FALSE
  note <- "kept_as_factor"

  if (is.numeric(x) || is.integer(x)) {
    use_numeric <- TRUE
    note <- "raw_numeric_or_integer"
  }

  if (force_numeric && parse_fraction >= parse_min_fraction) {
    use_numeric <- TRUE
    note <- "forced_numeric_parseable"
  }

  if (!force_factor && prefer_numeric && parse_fraction >= parse_min_fraction) {
    use_numeric <- TRUE
    note <- "preferred_numeric_parseable"
  }

  if (!force_factor && raw_n_unique > high_cardinality_threshold &&
      parse_fraction >= parse_min_fraction) {
    use_numeric <- TRUE
    note <- "high_cardinality_numeric_like"
  }

  if (force_factor) {
    use_numeric <- FALSE
    note <- "forced_factor"
  }

  if (use_numeric) {
    value <- as.numeric(parsed)
    encoding <- "numeric"
  } else {
    if (raw_n_unique > high_cardinality_threshold && collapse_high_cardinality) {
      value <- collapse_factor_levels(clean, max_levels = max_factor_levels)
      encoding <- "collapsed_factor"
      note <- paste0(
        "nonparseable_high_cardinality_factor_collapsed_to_top_",
        max_factor_levels,
        "_plus_other"
      )
    } else {
      value <- factor(clean)
      encoding <- "factor"
    }
  }

  final_n_unique <- length(unique(value[!is.na(value)]))
  final_class <- paste(class(value), collapse = "/")

  diagnostics <- data.frame(
    variable = variable,
    raw_class = raw_class,
    raw_n_unique = raw_n_unique,
    n_missing_raw = n_missing_raw,
    numeric_parse_fraction = as.numeric(parse_fraction),
    force_numeric = force_numeric,
    force_factor = force_factor,
    prefer_numeric = prefer_numeric,
    final_encoding = encoding,
    final_class = final_class,
    final_n_unique = final_n_unique,
    note = note,
    stringsAsFactors = FALSE
  )

  list(
    value = value,
    diagnostics = diagnostics
  )
}

rebuild_encoding_diagnostics_after_drop <- function(analysis, diagnostics) {

  diagnostics$n_missing_after_drop <- vapply(
    diagnostics$variable,
    function(v) {
      if (v %in% names(analysis)) {
        sum(is.na(analysis[[v]]))
      } else {
        NA_real_
      }
    },
    numeric(1)
  )

  diagnostics$final_n_unique_after_drop <- vapply(
    diagnostics$variable,
    function(v) {
      if (v %in% names(analysis)) {
        length(unique(analysis[[v]][!is.na(analysis[[v]])]))
      } else {
        NA_real_
      }
    },
    numeric(1)
  )

  diagnostics
}

prepare_data <- function(data, config = NULL) {

  enc_gender <- encode_ps_variable(data$gender, "gender", config, default_role = "factor")
  enc_age <- encode_ps_variable(data$age, "age", config, default_role = "numeric")
  enc_income <- encode_ps_variable(data$income, "income", config, default_role = "numeric")
  enc_income_missing <- encode_ps_variable(data$income.missing, "income.missing", config, default_role = "numeric")
  enc_race <- encode_ps_variable(data$race, "race", config, default_role = "factor")
  enc_education <- encode_ps_variable(data$education, "education", config, default_role = "factor")
  enc_smoking_ever <- encode_ps_variable(data$smoking.ever, "smoking.ever", config, default_role = "factor")
  enc_smoking_now <- encode_ps_variable(data$smoking.now, "smoking.now", config, default_role = "numeric")

  out <- data.frame(
    Z = ifelse(data$fish.level == "high", 1, 0),
    Y = log2(data$o.LBXTHG),
    gender = enc_gender$value,
    age = enc_age$value,
    income = enc_income$value,
    income.missing = enc_income_missing$value,
    race = enc_race$value,
    education = enc_education$value,
    smoking.ever = enc_smoking_ever$value,
    smoking.now = enc_smoking_now$value
  )

  diagnostics <- do.call(
    rbind,
    list(
      enc_gender$diagnostics,
      enc_age$diagnostics,
      enc_income$diagnostics,
      enc_income_missing$diagnostics,
      enc_race$diagnostics,
      enc_education$diagnostics,
      enc_smoking_ever$diagnostics,
      enc_smoking_now$diagnostics
    )
  )

  keep <- complete.cases(out)
  out <- out[keep, , drop = FALSE]

  diagnostics <- rebuild_encoding_diagnostics_after_drop(
    analysis = out,
    diagnostics = diagnostics
  )

  attr(out, "encoding_diagnostics") <- diagnostics

  out
}

add_ps_model_columns_to_encoding_diagnostics <- function(
    diagnostics,
    full_formula,
    data,
    config = NULL
) {

  if (is.null(diagnostics) || nrow(diagnostics) == 0) {
    return(data.frame())
  }

  full_cols <- colnames(model.matrix(full_formula, data))
  high_cardinality_threshold <- .config_value(config, "ps_high_cardinality_threshold", 15)

  diagnostics$model_n_columns <- NA_real_
  diagnostics$model_columns <- NA_character_
  diagnostics$model_encoding_status <- NA_character_

  for (i in seq_len(nrow(diagnostics))) {
    v <- diagnostics$variable[i]

    reduced_formula <- try(
      update(full_formula, paste(". ~ . -", v)),
      silent = TRUE
    )

    if (inherits(reduced_formula, "try-error")) {
      diagnostics$model_encoding_status[i] <- "term_not_in_formula_or_update_failed"
      next
    }

    reduced_cols <- colnames(model.matrix(reduced_formula, data))
    removed_cols <- setdiff(full_cols, reduced_cols)

    diagnostics$model_n_columns[i] <- length(removed_cols)
    diagnostics$model_columns[i] <- paste(removed_cols, collapse = " | ")

    diagnostics$model_encoding_status[i] <- ifelse(
      length(removed_cols) > high_cardinality_threshold,
      "high_dimensional_in_model",
      "low_dimensional_in_model"
    )
  }

  diagnostics
}
