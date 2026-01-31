# Benchmark: RcppICA vs fastICA Scaling Comparison
#
# This script benchmarks RcppICA against fastICA::fastICA at various
# matrix sizes to demonstrate the speed and memory advantages of
# Spectra's partial eigendecomposition.
#
# Usage: Rscript benchmarks/scaling_comparison.R

library(RcppICA)
library(fastICA)

cat("=== RcppICA vs fastICA Scaling Benchmark ===\n\n")

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Matrix sizes: (n_rows, n_cols, n_components)
# Keeping n/m ratio ~ 10:1
configs <- data.frame(
  n = c(5000, 10000, 20000, 30000, 40000, 50000),
  m = c(500,  1000,  2000,  3000,  4000,  5000),
  k = c(10,   15,    20,    25,    25,    30)
)

# Output paths
output_rds <- "benchmarks/scaling_results.rds"
output_png <- "benchmarks/scaling_comparison.png"
output_pdf <- "benchmarks/scaling_comparison.pdf"

# ---------------------------------------------------------------------------
# Helper: Generate non-Gaussian ICA data
# ---------------------------------------------------------------------------

generate_ica_data <- function(n, m, k, seed = 42) {
  set.seed(seed)

  # k truly independent non-Gaussian sources
  S <- matrix(0, n, k)
  for (i in seq_len(k)) {
    type <- i %% 5
    S[, i] <- switch(as.character(type),
      "0" = sin(seq(0, (2 + i) * pi, length.out = n)),
      "1" = sign(rnorm(n)) * rexp(n),
      "2" = runif(n) - 0.5,
      "3" = ((1:n %% (50 + i * 10)) - 25) / 25,
      "4" = rnorm(n)^3
    )
  }

  # Random mixing matrix
  A <- matrix(rnorm(m * k), m, k)

  # Mixed signals + tiny noise
  X <- S %*% t(A) + matrix(rnorm(n * m, sd = 0.01), n, m)

  list(X = X, S = S, A = A)
}

# ---------------------------------------------------------------------------
# Run Benchmarks
# ---------------------------------------------------------------------------

results <- data.frame()

for (i in seq_len(nrow(configs))) {
  cfg <- configs[i, ]
  cat(sprintf("Running %d x %d, k=%d ...\n", cfg$n, cfg$m, cfg$k))

  # Generate data
  data <- generate_ica_data(cfg$n, cfg$m, cfg$k)
  X <- data$X

  # RcppICA
  gc(reset = TRUE)
  t_rcpp <- system.time(r_rcpp <- RcppICA::fastICA(X, n.comp = cfg$k, seed = 42L))
  gc_rcpp <- gc()
  mem_rcpp <- max(gc_rcpp[, 6])  # max used (Vcells)
  rcpp_converged <- r_rcpp@misc$converged
  rcpp_iters <- r_rcpp@misc$iterations
  rm(r_rcpp); gc()

  # fastICA
  gc(reset = TRUE)
  t_fast <- system.time(r_fast <- fastICA::fastICA(X, n.comp = cfg$k, method = "C"))
  gc_fast <- gc()
  mem_fast <- max(gc_fast[, 6])  # max used (Vcells)
  rm(r_fast); gc()

  results <- rbind(results, data.frame(
    n = cfg$n,
    m = cfg$m,
    k = cfg$k,
    label = paste0(cfg$n / 1000, "Kx", cfg$m / 1000, "K"),
    elements_M = (cfg$n * cfg$m) / 1e6,
    time_rcpp = as.numeric(t_rcpp["elapsed"]),
    time_fast = as.numeric(t_fast["elapsed"]),
    mem_rcpp_vcells = mem_rcpp,
    mem_fast_vcells = mem_fast,
    converged = rcpp_converged,
    iters = rcpp_iters,
    stringsAsFactors = FALSE
  ))

  cat(sprintf("  RcppICA: %.2fs, %.1f GB | fastICA: %.2fs, %.1f GB | Speedup: %.1fx\n",
              t_rcpp["elapsed"], mem_rcpp * 8 / 1024^3,
              t_fast["elapsed"], mem_fast * 8 / 1024^3,
              t_fast["elapsed"] / max(t_rcpp["elapsed"], 0.001)))

  rm(X, data); gc()
}

# Convert Vcells to GB (each Vcell = 8 bytes)
results$mem_rcpp_gb <- results$mem_rcpp_vcells * 8 / 1024^3
results$mem_fast_gb <- results$mem_fast_vcells * 8 / 1024^3
results$speedup <- results$time_fast / results$time_rcpp
results$mem_ratio <- results$mem_fast_gb / results$mem_rcpp_gb

# ---------------------------------------------------------------------------
# Save Results
# ---------------------------------------------------------------------------

cat("\n--- Results Table ---\n")
print(results[, c("label", "k", "time_rcpp", "time_fast", "speedup",
                  "mem_rcpp_gb", "mem_fast_gb", "mem_ratio")])

saveRDS(results, output_rds)
cat(sprintf("\nResults saved to: %s\n", output_rds))

# ---------------------------------------------------------------------------
# Generate Plots
# ---------------------------------------------------------------------------

create_plot <- function(device_fn, path, width, height, ...) {
  device_fn(path, width = width, height = height, ...)
  par(mfrow = c(1, 2), mar = c(5, 4.5, 3, 1))

  # Panel 1: Execution Time
  y_max <- max(results$time_fast) * 1.1
  plot(results$elements_M, results$time_fast, type = "b", pch = 16, col = "firebrick",
       xlab = "Matrix Elements (millions)", ylab = "Time (seconds)",
       main = "Execution Time", ylim = c(0, y_max), lwd = 2, cex = 1.3)
  lines(results$elements_M, results$time_rcpp, type = "b", pch = 17, col = "steelblue",
        lwd = 2, cex = 1.3)
  legend("topleft", legend = c("fastICA", "RcppICA (Spectra)"),
         col = c("firebrick", "steelblue"), pch = c(16, 17), lwd = 2, bty = "n", cex = 1.1)

  for (j in seq_len(nrow(results))) {
    text(results$elements_M[j], results$time_rcpp[j],
         sprintf("%.1fx", results$speedup[j]), pos = 4, cex = 0.8, col = "steelblue")
  }

  # Panel 2: Peak Memory
  y_max_mem <- max(results$mem_fast_gb) * 1.1
  plot(results$elements_M, results$mem_fast_gb, type = "b", pch = 16, col = "firebrick",
       xlab = "Matrix Elements (millions)", ylab = "Peak Memory (GB)",
       main = "Peak Memory Usage", ylim = c(0, y_max_mem), lwd = 2, cex = 1.3)
  lines(results$elements_M, results$mem_rcpp_gb, type = "b", pch = 17, col = "steelblue",
        lwd = 2, cex = 1.3)
  legend("topleft", legend = c("fastICA", "RcppICA (Spectra)"),
         col = c("firebrick", "steelblue"), pch = c(16, 17), lwd = 2, bty = "n", cex = 1.1)

  for (j in seq_len(nrow(results))) {
    text(results$elements_M[j], results$mem_rcpp_gb[j],
         sprintf("%.1fx less", results$mem_ratio[j]), pos = 4, cex = 0.8, col = "steelblue")
  }

  dev.off()
}

# PNG
create_plot(png, output_png, width = 1200, height = 500, res = 120)
cat(sprintf("PNG saved to: %s\n", output_png))

# PDF
create_plot(pdf, output_pdf, width = 10, height = 5)
cat(sprintf("PDF saved to: %s\n", output_pdf))

cat("\nDone!\n")
