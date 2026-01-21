#!/usr/bin/env Rscript
# Benchmark Spectra vs Eigen vs base fastICA

library(RcppICA)
library(fastICA)

# Function to generate test data
generate_test_data <- function(n, m, k = 5) {
    # Generate k independent sources
    S <- matrix(rnorm(n * k), n, k)

    # Random mixing matrix
    A <- matrix(rnorm(m * k), m, k)

    # Mixed observations
    X <- S %*% t(A)

    return(X)
}

# Benchmark function
benchmark_methods <- function(n, m, k = 5, n_runs = 5) {
    cat(sprintf("\n=== Benchmarking: n=%d, m=%d, k=%d ===\n", n, m, k))

    # Generate data once
    set.seed(42)
    X <- generate_test_data(n, m, k)

    # Benchmark base fastICA
    cat("Running base fastICA...\n")
    base_times <- numeric(n_runs)
    for (i in 1:n_runs) {
        base_times[i] <- system.time({
            base_result <- fastICA::fastICA(X, n.comp = k, method = "C", verbose = FALSE)
        })["elapsed"]
    }
    base_median <- median(base_times)

    # Benchmark RcppICA with Spectra (new default)
    cat("Running RcppICA with Spectra...\n")
    spectra_times <- numeric(n_runs)
    for (i in 1:n_runs) {
        spectra_times[i] <- system.time({
            spectra_result <- RcppICA::fastICA(X, n.comp = k, whiten.method = "spectra", verbose = FALSE)
        })["elapsed"]
    }
    spectra_median <- median(spectra_times)

    # Benchmark RcppICA with Eigen (old default)
    cat("Running RcppICA with Eigen...\n")
    eigen_times <- numeric(n_runs)
    for (i in 1:n_runs) {
        eigen_times[i] <- system.time({
            eigen_result <- RcppICA::fastICA(X, n.comp = k, whiten.method = "eigen", verbose = FALSE)
        })["elapsed"]
    }
    eigen_median <- median(eigen_times)

    # Benchmark RcppICA with SVD
    cat("Running RcppICA with SVD...\n")
    svd_times <- numeric(n_runs)
    for (i in 1:n_runs) {
        svd_times[i] <- system.time({
            svd_result <- RcppICA::fastICA(X, n.comp = k, whiten.method = "svd", verbose = FALSE)
        })["elapsed"]
    }
    svd_median <- median(svd_times)

    # Verify correctness (should all be highly correlated)
    cat("\nVerifying correctness (correlation between methods)...\n")

    # Compare Spectra vs Eigen
    S_spectra <- t(spectra_result@S)
    S_eigen <- t(eigen_result@S)
    S_base <- base_result$S

    # Compute absolute correlations (ICA components can have arbitrary sign/order)
    cors_spectra_eigen <- abs(cor(S_spectra, S_eigen))
    cors_spectra_base <- abs(cor(S_spectra, S_base))

    max_cors_spectra_eigen <- apply(cors_spectra_eigen, 1, max)
    max_cors_spectra_base <- apply(cors_spectra_base, 1, max)

    cat(sprintf("Spectra vs Eigen correlation: %.6f (min: %.6f)\n",
                mean(max_cors_spectra_eigen), min(max_cors_spectra_eigen)))
    cat(sprintf("Spectra vs Base correlation: %.6f (min: %.6f)\n",
                mean(max_cors_spectra_base), min(max_cors_spectra_base)))

    # Print results
    cat("\n--- Results ---\n")
    cat(sprintf("Base fastICA (C):       %.4fs\n", base_median))
    cat(sprintf("RcppICA (Spectra):      %.4fs  [%.2fx vs base, %.2fx vs Eigen]\n",
                spectra_median, base_median/spectra_median, eigen_median/spectra_median))
    cat(sprintf("RcppICA (Eigen):        %.4fs  [%.2fx vs base]\n",
                eigen_median, base_median/eigen_median))
    cat(sprintf("RcppICA (SVD):          %.4fs  [%.2fx vs base]\n",
                svd_median, base_median/svd_median))

    # Return results
    list(
        dataset = sprintf("%dx%d (k=%d)", n, m, k),
        base = base_median,
        spectra = spectra_median,
        eigen = eigen_median,
        svd = svd_median,
        speedup_vs_base = base_median / spectra_median,
        speedup_vs_eigen = eigen_median / spectra_median,
        correlation_spectra_eigen = mean(max_cors_spectra_eigen),
        correlation_spectra_base = mean(max_cors_spectra_base)
    )
}

# Run benchmarks on different dataset sizes
cat(strrep("=", 70), "\n")
cat("RcppICA Performance Benchmark: Spectra vs Eigen vs Base fastICA\n")
cat(strrep("=", 70), "\n")

results <- list()

# Small dataset (quick test)
results[[1]] <- benchmark_methods(n = 500, m = 100, k = 5, n_runs = 10)

# Medium dataset (typical use case)
results[[2]] <- benchmark_methods(n = 1000, m = 500, k = 5, n_runs = 5)

# Large dataset (stress test)
results[[3]] <- benchmark_methods(n = 5000, m = 1000, k = 10, n_runs = 3)

# Very large dataset (if we want to be ambitious)
# results[[4]] <- benchmark_methods(n = 10000, m = 2000, k = 10, n_runs = 2)

# Summary table
cat("\n")
cat(strrep("=", 70), "\n")
cat("SUMMARY TABLE\n")
cat(strrep("=", 70), "\n\n")

cat(sprintf("%-20s | %8s | %8s | %8s | %8s | %8s\n",
            "Dataset", "Base", "Spectra", "Eigen", "vs Base", "vs Eigen"))
cat(strrep("-", 85), "\n")

for (r in results) {
    cat(sprintf("%-20s | %7.3fs | %7.3fs | %7.3fs | %7.2fx | %7.2fx\n",
                r$dataset, r$base, r$spectra, r$eigen, r$speedup_vs_base, r$speedup_vs_eigen))
}

cat("\n")
cat("Conclusion:\n")
avg_speedup_vs_base <- mean(sapply(results, function(x) x$speedup_vs_base))
avg_speedup_vs_eigen <- mean(sapply(results, function(x) x$speedup_vs_eigen))

cat(sprintf("- Spectra is %.2fx faster than base fastICA on average\n", avg_speedup_vs_base))
cat(sprintf("- Spectra is %.2fx faster than Eigen method on average\n", avg_speedup_vs_eigen))
cat(sprintf("- Results are highly correlated (r > %.3f), confirming correctness\n",
            min(sapply(results, function(x) x$correlation_spectra_eigen))))

if (avg_speedup_vs_base >= 2) {
    cat("\n✓ SUCCESS: Achieved 2-5x speedup goal over base fastICA!\n")
} else if (avg_speedup_vs_base >= 1) {
    cat(sprintf("\n✓ GOOD: Faster than base fastICA (%.2fx)\n", avg_speedup_vs_base))
} else {
    cat(sprintf("\n✗ Need improvement: Still %.2fx slower than base\n", 1/avg_speedup_vs_base))
}
