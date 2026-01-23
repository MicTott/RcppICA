#!/usr/bin/env Rscript
# Comprehensive Benchmark: RcppICA vs fastICA vs RcppML NMF
# Compares performance across different dataset sizes with statistical rigor

devtools::load_all()
library(fastICA)
library(RcppML)
library(Matrix)
library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(42)

# ===== Configuration =====
n_genes <- 5000
n_cells_sizes <- c(50000, 100000)
n_components <- 20
n_replicates <- 3
sparsity <- 0.10  # 10% non-zero (typical for scRNA-seq)

cat("=== Benchmark Configuration ===\n")
cat(sprintf("Genes: %d\n", n_genes))
cat(sprintf("Cell sizes: %s\n", paste(n_cells_sizes, collapse = ", ")))
cat(sprintf("Components: %d\n", n_components))
cat(sprintf("Replicates: %d\n", n_replicates))
cat(sprintf("Sparsity: %.1f%% non-zero\n\n", sparsity * 100))

# ===== Generate Simulated Data =====
generate_sparse_scrnaseq <- function(n_cells, n_genes, sparsity = 0.10) {
  # Simulate realistic scRNA-seq data
  # 1. Create sparse matrix with dropout
  nnz <- round(n_cells * n_genes * sparsity)

  i <- sample(1:n_cells, nnz, replace = TRUE)
  j <- sample(1:n_genes, nnz, replace = TRUE)

  # Log-normal-like distribution (typical for scRNA-seq)
  x <- rexp(nnz, rate = 0.1)

  X_sparse <- sparseMatrix(i = i, j = j, x = x, dims = c(n_cells, n_genes))

  # Also create dense version
  X_dense <- as.matrix(X_sparse)

  return(list(sparse = X_sparse, dense = X_dense))
}

# ===== Benchmark Functions =====

benchmark_method <- function(method_name, run_func, X, n_comp, n_rep = 3) {
  times <- numeric(n_rep)

  for (rep in 1:n_rep) {
    gc()  # Clean up memory before each run

    time <- system.time({
      result <- run_func(X, n_comp)
    })

    times[rep] <- time["elapsed"]
  }

  return(data.frame(
    method = method_name,
    replicate = 1:n_rep,
    time = times,
    mean_time = mean(times),
    sd_time = sd(times),
    se_time = sd(times) / sqrt(n_rep)
  ))
}

# ===== Method Implementations =====

# RcppICA - Dense
run_rcppica_dense <- function(X, k) {
  result <- fastICA(X, n.comp = k, alg.typ = "parallel",
                    fun = "logcosh", verbose = FALSE)
  return(result)
}

# RcppICA - Sparse
run_rcppica_sparse <- function(X, k) {
  result <- fastICA(X, n.comp = k, alg.typ = "parallel",
                    fun = "logcosh", verbose = FALSE)
  return(result)
}

# fastICA (R implementation)
run_fastica_r <- function(X, k) {
  result <- fastICA::fastICA(X, n.comp = k, alg.typ = "parallel",
                              fun = "logcosh", method = "R",
                              verbose = FALSE)
  return(result)
}

# fastICA (C implementation)
run_fastica_c <- function(X, k) {
  result <- fastICA::fastICA(X, n.comp = k, alg.typ = "parallel",
                              fun = "logcosh", method = "C",
                              verbose = FALSE)
  return(result)
}

# RcppML NMF (dense)
run_rcppml_nmf_dense <- function(X, k) {
  # NMF requires non-negative input
  X_nonneg <- X
  X_nonneg[X_nonneg < 0] <- 0

  result <- RcppML::nmf(X_nonneg, k = k,
                        tol = 1e-4, maxit = 100, verbose = FALSE)
  return(result)
}

# RcppML NMF (sparse)
run_rcppml_nmf_sparse <- function(X, k) {
  # NMF requires non-negative input
  X_nonneg <- X
  X_nonneg[X_nonneg < 0] <- 0
  
  result <- RcppML::nmf(X_nonneg, k = k,
                        tol = 1e-4, maxit = 100, verbose = FALSE)
  return(result)
}


# ===== Run Benchmarks =====

results_list <- list()
result_idx <- 1

cat("=== Running Benchmarks ===\n\n")

for (n_cells in n_cells_sizes) {
  cat(sprintf("Dataset size: %d cells × %d genes\n", n_cells, n_genes))

  # Generate data
  cat("  Generating data...")
  data <- generate_sparse_scrnaseq(n_cells, n_genes, sparsity)
  cat(" Done\n")

  # Memory footprint
  sparse_mb <- object.size(data$sparse) / 1024^2
  dense_mb <- object.size(data$dense) / 1024^2
  cat(sprintf("  Memory: Sparse = %.1f MB, Dense = %.1f MB (%.1fx)\n",
              sparse_mb, dense_mb, dense_mb / sparse_mb))

  # RcppICA - Dense
  cat("  Running RcppICA (dense)...")
  result <- benchmark_method("RcppICA_dense", run_rcppica_dense,
                             data$dense, n_components, n_replicates)
  result$n_cells <- n_cells
  result$n_genes <- n_genes
  results_list[[result_idx]] <- result
  result_idx <- result_idx + 1
  cat(sprintf(" %.3f ± %.3f sec\n", result$mean_time[1], result$sd_time[1]))

  # RcppICA - Sparse
  cat("  Running RcppICA (sparse)...")
  result <- benchmark_method("RcppICA_sparse", run_rcppica_sparse,
                             data$sparse, n_components, n_replicates)
  result$n_cells <- n_cells
  result$n_genes <- n_genes
  results_list[[result_idx]] <- result
  result_idx <- result_idx + 1
  cat(sprintf(" %.3f ± %.3f sec\n", result$mean_time[1], result$sd_time[1]))

  # fastICA (R)
  cat("  Running fastICA (R)...")
  result <- benchmark_method("fastICA_R", run_fastica_r,
                             data$dense, n_components, n_replicates)
  result$n_cells <- n_cells
  result$n_genes <- n_genes
  results_list[[result_idx]] <- result
  result_idx <- result_idx + 1
  cat(sprintf(" %.3f ± %.3f sec\n", result$mean_time[1], result$sd_time[1]))

  # fastICA (C)
  cat("  Running fastICA (C)...")
  result <- benchmark_method("fastICA_C", run_fastica_c,
                             data$dense, n_components, n_replicates)
  result$n_cells <- n_cells
  result$n_genes <- n_genes
  results_list[[result_idx]] <- result
  result_idx <- result_idx + 1
  cat(sprintf(" %.3f ± %.3f sec\n", result$mean_time[1], result$sd_time[1]))

  # RcppML NMF (dense)
  cat("  Running RcppML NMF (dense)...")
  result <- benchmark_method("RcppML_NMF", run_rcppml_nmf_dense,
                             data$dense, n_components, n_replicates)
  result$n_cells <- n_cells
  result$n_genes <- n_genes
  results_list[[result_idx]] <- result
  result_idx <- result_idx + 1
  cat(sprintf(" %.3f ± %.3f sec\n\n", result$mean_time[1], result$sd_time[1]))
  
  # RcppML NMF (sparse)
  cat("  Running RcppML NMF (dense)...")
  result <- benchmark_method("RcppML_NMF", run_rcppml_nmf_sparse,
                             data$sparse, n_components, n_replicates)
  result$n_cells <- n_cells
  result$n_genes <- n_genes
  results_list[[result_idx]] <- result
  result_idx <- result_idx + 1
  cat(sprintf(" %.3f ± %.3f sec\n\n", result$mean_time[1], result$sd_time[1]))

  gc()  # Clean up after each size
}

# ===== Combine Results =====
results_df <- bind_rows(results_list)

# Save raw results
write.csv(results_df, "benchmark_results.csv", row.names = FALSE)
cat("Raw results saved to: benchmark_results.csv\n\n")

# ===== Compute Summary Statistics =====
summary_stats <- results_df %>%
  group_by(method, n_cells, n_genes) %>%
  summarise(
    mean_time = mean(time),
    sd_time = sd(time),
    se_time = sd(time) / sqrt(n()),
    min_time = min(time),
    max_time = max(time),
    .groups = "drop"
  )

write.csv(summary_stats, "benchmark_summary.csv", row.names = FALSE)
cat("Summary statistics saved to: benchmark_summary.csv\n\n")

# ===== Print Summary Table =====
cat("=== Summary Results ===\n\n")
summary_wide <- summary_stats %>%
  select(method, n_cells, mean_time, se_time) %>%
  mutate(time_str = sprintf("%.3f ± %.3f", mean_time, se_time)) %>%
  select(method, n_cells, time_str) %>%
  pivot_wider(names_from = n_cells, values_from = time_str,
              names_prefix = "cells_")

print(summary_wide, n = 100)

# ===== Compute Speedup Relative to fastICA_C =====
cat("\n=== Speedup vs fastICA (C) ===\n\n")

baseline <- summary_stats %>%
  filter(method == "fastICA_C") %>%
  select(n_cells, baseline_time = mean_time)

speedup_df <- summary_stats %>%
  left_join(baseline, by = "n_cells") %>%
  mutate(speedup = baseline_time / mean_time) %>%
  select(method, n_cells, speedup) %>%
  pivot_wider(names_from = n_cells, values_from = speedup,
              names_prefix = "cells_")

print(speedup_df, n = 100)

# ===== Generate Plots =====
cat("\n=== Generating Plots ===\n")

# Color palette
method_colors <- c(
  "RcppICA_dense" = "#E64B35FF",
  "RcppICA_sparse" = "#00A087FF",
  "fastICA_R" = "#3C5488FF",
  "fastICA_C" = "#F39B7FFF",
  "RcppML_NMF" = "#8491B4FF"
)

# Method labels
method_labels <- c(
  "RcppICA_dense" = "RcppICA (dense)",
  "RcppICA_sparse" = "RcppICA (sparse)",
  "fastICA_R" = "fastICA (R)",
  "fastICA_C" = "fastICA (C)",
  "RcppML_NMF" = "RcppML NMF"
)

# Plot 1: Runtime vs Dataset Size (with error bars)
p1 <- ggplot(summary_stats, aes(x = n_cells, y = mean_time,
                                 color = method, group = method)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_time - se_time, ymax = mean_time + se_time),
                width = 0.1 * max(n_cells_sizes), alpha = 0.7) +
  scale_color_manual(values = method_colors, labels = method_labels) +
  scale_x_continuous(trans = "log10",
                     breaks = n_cells_sizes,
                     labels = scales::comma) +
  scale_y_continuous(trans = "log10") +
  labs(title = "Runtime Comparison: ICA and NMF Methods",
       subtitle = sprintf("%d genes × varying cells, %d components, %d replicates",
                          n_genes, n_components, n_replicates),
       x = "Number of Cells (log scale)",
       y = "Runtime (seconds, log scale)",
       color = "Method") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        panel.grid.minor = element_blank())

ggsave("benchmark_runtime_log.png", p1, width = 10, height = 8, dpi = 300)
cat("  Saved: benchmark_runtime_log.png\n")

# Plot 2: Linear scale (for smaller datasets)
summary_stats_small <- summary_stats %>% filter(n_cells <= 10000)

p2 <- ggplot(summary_stats_small, aes(x = n_cells, y = mean_time,
                                       color = method, group = method)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_time - se_time, ymax = mean_time + se_time),
                width = 200, alpha = 0.7) +
  geom_ribbon(aes(ymin = mean_time - se_time, ymax = mean_time + se_time,
                  fill = method), alpha = 0.2, color = NA) +
  scale_color_manual(values = method_colors, labels = method_labels) +
  scale_fill_manual(values = method_colors, labels = method_labels) +
  scale_x_continuous(breaks = c(1000, 5000, 10000),
                     labels = scales::comma) +
  labs(title = "Runtime Comparison (Linear Scale)",
       subtitle = sprintf("Small to medium datasets (%d genes, %d components)",
                          n_genes, n_components),
       x = "Number of Cells",
       y = "Runtime (seconds)",
       color = "Method",
       fill = "Method") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        panel.grid.minor = element_blank())

ggsave("benchmark_runtime_linear.png", p2, width = 10, height = 8, dpi = 300)
cat("  Saved: benchmark_runtime_linear.png\n")

# Plot 3: Speedup relative to fastICA (C)
speedup_plot_df <- summary_stats %>%
  left_join(baseline, by = "n_cells") %>%
  mutate(speedup = baseline_time / mean_time,
         # Propagate error for speedup (simplified)
         speedup_se = speedup * sqrt((se_time / mean_time)^2))

p3 <- ggplot(speedup_plot_df, aes(x = n_cells, y = speedup,
                                   color = method, group = method)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", size = 1) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = speedup - speedup_se, ymax = speedup + speedup_se),
                width = 0.1 * max(n_cells_sizes), alpha = 0.7) +
  scale_color_manual(values = method_colors, labels = method_labels) +
  scale_x_continuous(trans = "log10",
                     breaks = n_cells_sizes,
                     labels = scales::comma) +
  labs(title = "Speedup Relative to fastICA (C)",
       subtitle = "Values > 1 indicate faster than baseline",
       x = "Number of Cells (log scale)",
       y = "Speedup Factor",
       color = "Method") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        panel.grid.minor = element_blank())

ggsave("benchmark_speedup.png", p3, width = 10, height = 8, dpi = 300)
cat("  Saved: benchmark_speedup.png\n")

# Plot 4: Individual replicate points (to show variability)
p4 <- ggplot(results_df, aes(x = n_cells, y = time, color = method)) +
  geom_point(alpha = 0.5, size = 2, position = position_jitter(width = 0.02)) +
  geom_line(data = summary_stats, aes(y = mean_time, group = method),
            size = 1.2) +
  scale_color_manual(values = method_colors, labels = method_labels) +
  scale_x_continuous(trans = "log10",
                     breaks = n_cells_sizes,
                     labels = scales::comma) +
  scale_y_continuous(trans = "log10") +
  labs(title = "Runtime Variability Across Replicates",
       subtitle = sprintf("%d replicates per method/size combination",
                          n_replicates),
       x = "Number of Cells (log scale)",
       y = "Runtime (seconds, log scale)",
       color = "Method") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        legend.direction = "vertical")

ggsave("benchmark_variability.png", p4, width = 10, height = 8, dpi = 300)
cat("  Saved: benchmark_variability.png\n")

# Plot 5: Bar chart for specific size (e.g., 10K cells)
bar_data <- summary_stats %>%
  filter(n_cells == 10000) %>%
  arrange(mean_time)

p5 <- ggplot(bar_data, aes(x = reorder(method, mean_time), y = mean_time,
                           fill = method)) +
  geom_col(alpha = 0.8) +
  geom_errorbar(aes(ymin = mean_time - se_time, ymax = mean_time + se_time),
                width = 0.3) +
  scale_fill_manual(values = method_colors, labels = method_labels) +
  scale_x_discrete(labels = method_labels) +
  labs(title = "Runtime Comparison at 10,000 Cells",
       subtitle = sprintf("%d genes, %d components", n_genes, n_components),
       x = "",
       y = "Runtime (seconds)") +
  coord_flip() +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank())

ggsave("benchmark_10k_cells.png", p5, width = 10, height = 6, dpi = 300)
cat("  Saved: benchmark_10k_cells.png\n")

# ===== Memory Comparison Plot =====
cat("\n=== Generating Memory Comparison ===\n")

memory_data <- data.frame(
  n_cells = n_cells_sizes,
  sparse_mb = sapply(n_cells_sizes, function(n) {
    nnz <- round(n * n_genes * sparsity)
    # dgCMatrix: 3 vectors (i, p, x) ≈ 3 * nnz * 8 bytes
    (nnz * 3 * 8) / 1024^2
  }),
  dense_mb = sapply(n_cells_sizes, function(n) {
    (n * n_genes * 8) / 1024^2
  })
)

memory_long <- memory_data %>%
  pivot_longer(cols = c(sparse_mb, dense_mb),
               names_to = "type", values_to = "memory_mb") %>%
  mutate(type = factor(type, levels = c("dense_mb", "sparse_mb"),
                       labels = c("Dense Matrix", "Sparse Matrix")))

p6 <- ggplot(memory_long, aes(x = n_cells, y = memory_mb,
                               color = type, group = type)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_x_continuous(trans = "log10",
                     breaks = n_cells_sizes,
                     labels = scales::comma) +
  scale_y_continuous(trans = "log10",
                     labels = scales::comma) +
  scale_color_manual(values = c("Dense Matrix" = "#E64B35FF",
                                 "Sparse Matrix" = "#00A087FF")) +
  labs(title = "Memory Footprint: Sparse vs Dense",
       subtitle = sprintf("%d genes, %.0f%% sparsity",
                          n_genes, (1 - sparsity) * 100),
       x = "Number of Cells (log scale)",
       y = "Memory (MB, log scale)",
       color = "Matrix Type") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

ggsave("benchmark_memory.png", p6, width = 10, height = 6, dpi = 300)
cat("  Saved: benchmark_memory.png\n")

# ===== Final Summary =====
cat("\n=== Benchmark Complete ===\n\n")
cat("Files generated:\n")
cat("  - benchmark_results.csv (raw data)\n")
cat("  - benchmark_summary.csv (summary statistics)\n")
cat("  - benchmark_runtime_log.png (log-log plot)\n")
cat("  - benchmark_runtime_linear.png (linear plot, small datasets)\n")
cat("  - benchmark_speedup.png (speedup vs fastICA C)\n")
cat("  - benchmark_variability.png (replicate variability)\n")
cat("  - benchmark_10k_cells.png (bar chart at 10K cells)\n")
cat("  - benchmark_memory.png (memory comparison)\n\n")

# Print key findings
cat("=== Key Findings ===\n\n")

fastest_overall <- summary_stats %>%
  group_by(n_cells) %>%
  slice_min(mean_time, n = 1) %>%
  select(n_cells, method, mean_time)

cat("Fastest method per dataset size:\n")
print(fastest_overall)

cat("\n")

# Sparse vs Dense speedup
sparse_dense_comparison <- summary_stats %>%
  filter(method %in% c("RcppICA_sparse", "RcppICA_dense")) %>%
  select(method, n_cells, mean_time) %>%
  pivot_wider(names_from = method, values_from = mean_time) %>%
  mutate(speedup = RcppICA_dense / RcppICA_sparse)

cat("RcppICA Sparse vs Dense speedup:\n")
print(sparse_dense_comparison)

cat("\nDone!\n")
