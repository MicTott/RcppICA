#!/usr/bin/env Rscript
# Profile to understand where time is spent

library(RcppICA)

# Generate large test data
set.seed(42)
n <- 5000
m <- 1000
k <- 10

S <- matrix(rnorm(n * k), n, k)
A <- matrix(rnorm(m * k), m, k)
X <- S %*% t(A)

cat("Profiling Spectra whitening on 5000x1000 data (k=10)\n")
cat(strrep("=", 50), "\n\n")

# Test just whitening performance
cat("Testing whitening methods only (no ICA)...\n\n")

# We'll need to access the C++ whitening functions directly
# For now, let's time full ICA but with very few iterations
# to isolate whitening time

cat("Spectra method (maxit=1 to isolate whitening):\n")
t1 <- system.time({
    r1 <- fastICA(X, n.comp = k, whiten.method = "spectra", maxit = 1, verbose = FALSE)
})
print(t1)

cat("\nEigen method (maxit=1 to isolate whitening):\n")
t2 <- system.time({
    r2 <- fastICA(X, n.comp = k, whiten.method = "eigen", maxit = 1, verbose = FALSE)
})
print(t2)

cat("\nSVD method (maxit=1 to isolate whitening):\n")
t3 <- system.time({
    r3 <- fastICA(X, n.comp = k, whiten.method = "svd", maxit = 1, verbose = FALSE)
})
print(t3)

cat("\n\nFull ICA run:\n")
cat("Spectra (full ICA):\n")
t4 <- system.time({
    r4 <- fastICA(X, n.comp = k, whiten.method = "spectra", verbose = FALSE)
})
print(t4)

cat(sprintf("\nWhitening overhead: %.3fs\n", t1["elapsed"]))
cat(sprintf("ICA iterations time: %.3fs\n", t4["elapsed"] - t1["elapsed"]))

# Compare with base fastICA
cat("\nBase fastICA (full ICA):\n")
t5 <- system.time({
    r5 <- fastICA::fastICA(X, n.comp = k, verbose = FALSE, method = "C")
})
print(t5)

cat("\n\nAnalysis:\n")
cat(sprintf("Spectra whitening only: %.3fs\n", t1["elapsed"]))
cat(sprintf("Eigen whitening only: %.3fs\n", t2["elapsed"]))
cat(sprintf("SVD whitening only: %.3fs\n", t3["elapsed"]))
cat(sprintf("Spectra speedup vs Eigen: %.2fx\n", t2["elapsed"] / t1["elapsed"]))
cat(sprintf("Spectra speedup vs SVD: %.2fx\n", t3["elapsed"] / t1["elapsed"]))
