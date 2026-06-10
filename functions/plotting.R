# ============================================================
# Plotting functions
# ============================================================

library(ggplot2)

plot_bootstrap_curve <- function(results, R2_star = NA_real_) {

  p <- ggplot(results, aes(x = R2)) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      fill = "grey60",
      alpha = 0.5
    ) +
    geom_line(aes(y = lower)) +
    geom_line(aes(y = upper)) +
    geom_hline(
      yintercept = 0,
      color = "red",
      linetype = "dashed"
    ) +
    facet_wrap(~stage) +
    theme_bw() +
    labs(
      title = "Variance-based bootstrap sensitivity inference",
      x = expression(R^2),
      y = "Bootstrap CI"
    )

  if (is.finite(R2_star)) {
    p <- p +
      geom_vline(
        xintercept = R2_star,
        linetype = "dotted",
        color = "red"
      )
  }

  p
}

normalize_benchmark_method_label <- function(x) {

  x <- as.character(x)

  out <- ifelse(
    grepl("Qbal", x, ignore.case = TRUE),
    "MSM (Qbal)",
    ifelse(
      grepl("Corr", x, ignore.case = TRUE),
      "VBM, w/ Corr.",
      ifelse(
        grepl("MSM", x, ignore.case = TRUE),
        "MSM",
        ifelse(
          grepl("VBM", x, ignore.case = TRUE),
          "VBM",
          x
        )
      )
    )
  )

  factor(
    out,
    levels = c("MSM", "MSM (Qbal)", "VBM", "VBM, w/ Corr.")
  )
}

plot_covariate_benchmark <- function(
    df,
    tau_hat = NA_real_,
    config = NULL
) {

  plot_df <- df
  plot_df$method_label <- normalize_benchmark_method_label(plot_df$method)

  if (!is.null(config$benchmark_plot_labels)) {
    label_map <- config$benchmark_plot_labels
    plot_df$variable_label <- ifelse(
      plot_df$variable %in% names(label_map),
      unname(label_map[plot_df$variable]),
      plot_df$variable
    )

    ordered_vars <- names(label_map)[names(label_map) %in% unique(plot_df$variable)]
    ordered_labels <- unname(label_map[ordered_vars])

    plot_df$variable_label <- factor(
      plot_df$variable_label,
      levels = ordered_labels
    )
  } else {
    plot_df$variable_label <- factor(
      plot_df$variable,
      levels = unique(plot_df$variable)
    )
  }

  p <- ggplot(
    plot_df,
    aes(
      x = variable_label,
      ymin = lower,
      ymax = upper,
      color = method_label
    )
  ) +
    geom_linerange(
      position = position_dodge(width = 0.75),
      linewidth = 1.05
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1)
    ) +
    labs(
      title = "Covariate benchmark: four sensitivity-model intervals",
      subtitle = "Intervals use shared percentile-bootstrap pools at fixed original-sample benchmark parameters.",
      x = "Benchmark covariate",
      y = "ATT interval",
      color = "Model"
    )

  if (is.finite(tau_hat)) {
    p <- p +
      geom_hline(
        yintercept = tau_hat,
        linetype = "dotted"
      )
  }

  p
}
