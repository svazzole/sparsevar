#!/usr/bin/env Rscript
#
# Simulation Comparison Script
# Compares different VAR estimation methods (ENET, SCAD, MCP) with and without
# thresholding across multiple simulations.
#

# Load the sparsevar package
library(sparsevar)
library(ggplot2)
library(tidyr)
library(dplyr)

# Set seed for reproducibility
set.seed(42)

# Simulation parameters
n_simulations <- 100        # Number of simulations
n_vars <- 30               # Dimension of VAR
p <- 1                     # VAR order (lag)
n_obs <- 250               # Number of observations
sparsity <- 0.05           # Sparsity of the VAR matrix
rho <- 0.5                 # Base covariance value

# Define estimation methods to compare
methods <- c(
  "ENET",
  "ENET_threshold",
  "SCAD",
  "SCAD_threshold",
  "MCP",
  "MCP_threshold"
)

# Initialize results data frame
results <- data.frame(
  simulation = integer(),
  method = character(),
  accuracy = numeric(),
  tpr = numeric(),
  tnr = numeric(),
  frobenius = numeric(),
  stringsAsFactors = FALSE
)

# Progress message
cat("Starting simulation comparison\n")
cat(sprintf("Parameters: n=%d, p=%d, nobs=%d, sparsity=%.2f, n_sims=%d\n",
            n_vars, p, n_obs, sparsity, n_simulations))
cat(paste(rep("-", 60), collapse = ""), "\n")

# Run simulations
for (sim in seq_len(n_simulations)) {
  cat(sprintf("Simulation %d/%d\n", sim, n_simulations))

  # Generate VAR data
  sim_data <- tryCatch({
    simulate_var(
      n = n_vars,
      p = p,
      nobs = n_obs,
      rho = rho,
      sparsity = sparsity,
      method = "normal",
      covariance = "Toeplitz"
    )
  }, error = function(e) {
    cat(sprintf("  Error in simulation: %s\n", e$message))
    return(NULL)
  })

  if (is.null(sim_data)) next

  # True coefficient matrix (first lag)
  true_matrix <- sim_data$A[[1]]

  # Fit each method
  for (method in methods) {
    # Parse method name
    use_threshold <- grepl("_threshold", method)
    penalty <- gsub("_threshold", "", method)

    cat(sprintf("  Fitting %s (threshold=%s)... ", penalty, use_threshold))

    fit_result <- tryCatch({
      fit_var(
        data = sim_data$series,
        p = p,
        penalty = penalty,
        method = "cv",
        threshold = use_threshold,
        nfolds = 5  # Fewer folds to speed up
      )
    }, error = function(e) {
      cat(sprintf("Error: %s\n", e$message))
      return(NULL)
    })

    if (is.null(fit_result)) next

    # Estimated coefficient matrix
    estimated_matrix <- fit_result$A[[1]]

    # Compute metrics
    acc <- accuracy(true_matrix, estimated_matrix)
    tpr <- true_positive_rate(true_matrix, estimated_matrix)
    tnr <- true_negative_rate(true_matrix, estimated_matrix)
    frob <- frob_norm(true_matrix - estimated_matrix)

    cat(sprintf("ACC=%.3f, TPR=%.3f, TNR=%.3f, Frob=%.3f\n",
                acc, tpr, tnr, frob))

    # Store results
    results <- rbind(results, data.frame(
      simulation = sim,
      method = method,
      accuracy = acc,
      tpr = tpr,
      tnr = tnr,
      frobenius = frob,
      stringsAsFactors = FALSE
    ))
  }
}

cat(paste(rep("-", 60), collapse = ""), "\n")
cat("Simulations completed. Creating plots...\n")

# Convert method to factor with proper ordering
results$method <- factor(results$method, levels = methods)

# Reshape data for plotting
results_long <- results %>%
  pivot_longer(
    cols = c(accuracy, tpr, tnr, frobenius),
    names_to = "metric",
    values_to = "value"
  )

# Create nice metric labels
metric_labels <- c(
  "accuracy" = "Accuracy",
  "tpr" = "True Positive Rate",
  "tnr" = "True Negative Rate",
  "frobenius" = "Frobenius Norm"
)

results_long$metric <- factor(
  results_long$metric,
  levels = names(metric_labels),
  labels = metric_labels
)

# Create method labels for better readability
method_labels <- c(
  "ENET" = "ENET",
  "ENET_threshold" = "ENET\n(thresh)",
  "SCAD" = "SCAD",
  "SCAD_threshold" = "SCAD\n(thresh)",
  "MCP" = "MCP",
  "MCP_threshold" = "MCP\n(thresh)"
)

# Create the combined boxplot
p_combined <- ggplot(results_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.5) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  scale_x_discrete(labels = method_labels) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Comparison of VAR Estimation Methods",
    subtitle = sprintf("VAR(%d), n=%d, T=%d, sparsity=%.0f%%, %d simulations",
                       p, n_vars, n_obs, sparsity * 100, n_simulations),
    x = "Estimation Method",
    y = "Value"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 8),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Save the combined plot
ggsave(
  "simulation_comparison_combined.pdf",
  plot = p_combined,
  width = 10,
  height = 8
)
cat("Saved: simulation_comparison_combined.pdf\n")

# Create individual metric plots for more detail
create_metric_plot <- function(data, metric_name, title, y_label) {
  metric_data <- data[data$metric == metric_name, ]

  ggplot(metric_data, aes(x = method, y = value, fill = method)) +
    geom_boxplot(outlier.shape = 21, outlier.size = 2) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    scale_x_discrete(labels = method_labels) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = title,
      subtitle = sprintf("VAR(%d), n=%d, T=%d, sparsity=%.0f%%",
                         p, n_vars, n_obs, sparsity * 100),
      x = "Estimation Method",
      y = y_label
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
}

# Create and save individual plots
p_acc <- create_metric_plot(results_long, "Accuracy",
                            "Accuracy Comparison", "Accuracy")
ggsave("simulation_accuracy.pdf", plot = p_acc, width = 8, height = 6)
cat("Saved: simulation_accuracy.pdf\n")

p_tpr <- create_metric_plot(results_long, "True Positive Rate",
                            "True Positive Rate Comparison",
                            "TPR (Sensitivity)")
ggsave("simulation_tpr.pdf", plot = p_tpr, width = 8, height = 6)
cat("Saved: simulation_tpr.pdf\n")

p_tnr <- create_metric_plot(results_long, "True Negative Rate",
                            "True Negative Rate Comparison",
                            "TNR (Specificity)")
ggsave("simulation_tnr.pdf", plot = p_tnr, width = 8, height = 6)
cat("Saved: simulation_tnr.pdf\n")

p_frob <- create_metric_plot(results_long, "Frobenius Norm",
                             "Frobenius Norm Comparison",
                             "Frobenius Norm (lower is better)")
ggsave("simulation_frobenius.pdf", plot = p_frob, width = 8, height = 6)
cat("Saved: simulation_frobenius.pdf\n")

# Print summary statistics
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("SUMMARY STATISTICS\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

summary_stats <- results %>%
  group_by(method) %>%
  summarise(
    n = n(),
    acc_mean = mean(accuracy, na.rm = TRUE),
    acc_sd = sd(accuracy, na.rm = TRUE),
    tpr_mean = mean(tpr, na.rm = TRUE),
    tpr_sd = sd(tpr, na.rm = TRUE),
    tnr_mean = mean(tnr, na.rm = TRUE),
    tnr_sd = sd(tnr, na.rm = TRUE),
    frob_mean = mean(frobenius, na.rm = TRUE),
    frob_sd = sd(frobenius, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats, width = Inf)

# Save results to CSV
write.csv(results, "simulation_results.csv", row.names = FALSE)
cat("\nSaved: simulation_results.csv\n")

write.csv(summary_stats, "simulation_summary.csv", row.names = FALSE)
cat("Saved: simulation_summary.csv\n")

cat("\nDone!\n")
