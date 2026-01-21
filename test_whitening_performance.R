# Test whitening performance - this is the bottleneck

library(RcppICA)

set.seed(42)
n <- 1000
m <- 500
X <- matrix(rnorm(n * m), n, m)

cat("Matrix dimensions:", n, "x", m, "\n")
cat("Testing whitening performance...\n\n")

# Test SVD whitening (our current approach)
cat("SVD whitening:\n")
time_svd <- system.time({
  for (i in 1:10) {
    result_svd <- RcppICA:::cpp_whiten(X, 5L, 0L)  # SVD method
  }
})
print(time_svd)

# Test Eigen whitening
cat("\nEigen whitening:\n")
time_eigen <- system.time({
  for (i in 1:10) {
    result_eigen <- RcppICA:::cpp_whiten(X, 5L, 1L)  # Eigen method
  }
})
print(time_eigen)

# Test what base R SVD does
cat("\nBase R svd():\n")
time_r_svd <- system.time({
  for (i in 1:10) {
    X_centered <- scale(X, center = TRUE, scale = FALSE)
    s <- svd(X_centered, nu = 5, nv = 5)  # Only compute first 5 components
  }
})
print(time_r_svd)

# Test what base R eigen does
cat("\nBase R eigen() on covariance:\n")
time_r_eigen <- system.time({
  for (i in 1:10) {
    X_centered <- scale(X, center = TRUE, scale = FALSE)
    C <- crossprod(X_centered) / (n - 1)
    e <- eigen(C, symmetric = TRUE)
  }
})
print(time_r_eigen)

cat("\nSummary:\n")
cat("SVD whitening:    ", round(time_svd["elapsed"]/10, 4), "s per call\n")
cat("Eigen whitening:  ", round(time_eigen["elapsed"]/10, 4), "s per call\n")
cat("Base R svd():     ", round(time_r_svd["elapsed"]/10, 4), "s per call\n")
cat("Base R eigen():   ", round(time_r_eigen["elapsed"]/10, 4), "s per call\n")
