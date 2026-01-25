#!/usr/bin/env Rscript
#
# ENET Threshold Comparison Script
# Compares ENET estimation with and without thresholding across multiple
# simulations.
#

# Load the sparsevar package
library(sparsevar)
library(ggplot2)
library(tidyr)
library(dplyr)

# Set seed for reproducibility
set.seed(42)

# Simulation parameters
n_simulations <- 100       # Number of simulations
n_vars <- 10               # Dimension of VAR
p <- 1                     # VAR order (lag)
n_obs <- 200               # Number of observations
sparsity <- 0.8            # Sparsity of the VAR matrix
rho <- 0.5                 # Base covariance value

# Define estimation methods to compare
methods <- c(
  "ENET",
  "ENET_threshold"
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
cat("Starting ENET threshold comparison\n")
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

    cat(sprintf("  Fitting ENET (threshold=%s)... ", use_threshold))

    fit_result <- tryCatch({
      fit_var(
        data = sim_data$series,
        p = p,
        penalty = "ENET",
        method = "cv",
        threshold = use_threshold,
        nfolds = 5
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
  "ENET_threshold" = "ENET\n(threshold)"
)

# Color palette for two methods
method_colors <- c("ENET" = "#66c2a5", "ENET_threshold" = "#fc8d62")

# Create the combined boxplot
p_combined <- ggplot(results_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.5) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  scale_x_discrete(labels = method_labels) +
  scale_fill_manual(values = method_colors) +
  labs(
    title = "ENET: Effect of Thresholding on Estimation Quality",
    subtitle = sprintf("VAR(%d), n=%d, T=%d, sparsity=%.0f%%, %d simulations",
                       p, n_vars, n_obs, sparsity * 100, n_simulations),
    x = "Estimation Method",
    y = "Value"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Save the combined plot
ggsave(
  "enet_comparison_combined.pdf",
  plot = p_combined,
  width = 9,
  height = 7
)
cat("Saved: enet_comparison_combined.pdf\n")

# Create individual metric plots for more detail
create_metric_plot <- function(data, metric_name, title, y_label) {
  metric_data <- data[data$metric == metric_name, ]

  ggplot(metric_data, aes(x = method, y = value, fill = method)) +
    geom_boxplot(outlier.shape = 21, outlier.size = 2, width = 0.5) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
    scale_x_discrete(labels = method_labels) +
    scale_fill_manual(values = method_colors) +
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
      axis.text.x = element_text(size = 11),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
}

# Create and save individual plots
p_acc <- create_metric_plot(results_long, "Accuracy",
                            "ENET Accuracy: With vs Without Threshold",
                            "Accuracy")
ggsave("enet_accuracy.pdf", plot = p_acc, width = 6, height = 5)
cat("Saved: enet_accuracy.pdf\n")

p_tpr <- create_metric_plot(results_long, "True Positive Rate",
                            "ENET True Positive Rate: With vs Without Threshold",
                            "TPR (Sensitivity)")
ggsave("enet_tpr.pdf", plot = p_tpr, width = 6, height = 5)
cat("Saved: enet_tpr.pdf\n")

p_tnr <- create_metric_plot(results_long, "True Negative Rate",
                            "ENET True Negative Rate: With vs Without Threshold",
                            "TNR (Specificity)")
ggsave("enet_tnr.pdf", plot = p_tnr, width = 6, height = 5)
cat("Saved: enet_tnr.pdf\n")

p_frob <- create_metric_plot(results_long, "Frobenius Norm",
                             "ENET Frobenius Norm: With vs Without Threshold",
                             "Frobenius Norm (lower is better)")
ggsave("enet_frobenius.pdf", plot = p_frob, width = 6, height = 5)
cat("Saved: enet_frobenius.pdf\n")

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

# Perform statistical tests (paired t-tests since same simulation data)
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("STATISTICAL TESTS (Paired t-tests)\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# Reshape for paired comparison
enet_results <- results[results$method == "ENET", ]
enet_thresh_results <- results[results$method == "ENET_threshold", ]

# Merge by simulation
paired_data <- merge(
  enet_results[, c("simulation", "accuracy", "tpr", "tnr", "frobenius")],
  enet_thresh_results[, c("simulation", "accuracy", "tpr", "tnr", "frobenius")],
  by = "simulation",
  suffixes = c("_enet", "_thresh")
)

if (nrow(paired_data) > 1) {
  cat("Accuracy:\n")
  t_acc <- t.test(paired_data$accuracy_enet, paired_data$accuracy_thresh,
                  paired = TRUE)
  cat(sprintf("  ENET mean: %.4f, ENET+thresh mean: %.4f, p-value: %.4f\n",
              mean(paired_data$accuracy_enet),
              mean(paired_data$accuracy_thresh),
              t_acc$p.value))

  cat("\nTrue Positive Rate:\n")
  t_tpr <- t.test(paired_data$tpr_enet, paired_data$tpr_thresh, paired = TRUE)
  cat(sprintf("  ENET mean: %.4f, ENET+thresh mean: %.4f, p-value: %.4f\n",
              mean(paired_data$tpr_enet),
              mean(paired_data$tpr_thresh),
              t_tpr$p.value))

  cat("\nTrue Negative Rate:\n")
  t_tnr <- t.test(paired_data$tnr_enet, paired_data$tnr_thresh, paired = TRUE)
  cat(sprintf("  ENET mean: %.4f, ENET+thresh mean: %.4f, p-value: %.4f\n",
              mean(paired_data$tnr_enet),
              mean(paired_data$tnr_thresh),
              t_tnr$p.value))

  cat("\nFrobenius Norm:\n")
  t_frob <- t.test(paired_data$frobenius_enet, paired_data$frobenius_thresh,
                   paired = TRUE)
  cat(sprintf("  ENET mean: %.4f, ENET+thresh mean: %.4f, p-value: %.4f\n",
              mean(paired_data$frobenius_enet),
              mean(paired_data$frobenius_thresh),
              t_frob$p.value))
}

# Save results to CSV
write.csv(results, "enet_comparison_results.csv", row.names = FALSE)
cat("\nSaved: enet_comparison_results.csv\n")

write.csv(summary_stats, "enet_comparison_summary.csv", row.names = FALSE)
cat("Saved: enet_comparison_summary.csv\n")

cat("\nDone!\n")
