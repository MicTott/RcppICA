# Performance profiling script to identify bottlenecks in RcppICA

library(RcppICA)

set.seed(12345)

# Create test data
n <- 1000
n_genes <- 500
n_comp <- 5

# Simulate data
S_true <- matrix(0, n, n_comp)
for (i in 1:n_comp) {
  S_true[, i] <- sin(seq(0, 10*pi, length.out = n)) + rnorm(n, 0, 0.1)
}
A_true <- matrix(runif(n_genes * n_comp, 0.3, 0.7), n_genes, n_comp)
X_expression <- t(S_true %*% t(A_true))

cat("Data dimensions: ", dim(X_expression), "\n")

# Time RcppICA
cat("\n=== Timing RcppICA ===\n")
time_rcpp <- system.time({
  result_rcpp <- fastICA(X_expression, n.comp = n_comp, maxit = 100,
                          alg.typ = "parallel", n.threads = 1, seed = 42)
})
print(time_rcpp)

# Time base fastICA if available
if (requireNamespace("fastICA", quietly = TRUE)) {
  cat("\n=== Timing base fastICA ===\n")
  time_base <- system.time({
    result_base <- fastICA::fastICA(X_expression, n.comp = n_comp, maxit = 100,
                                     fun = "logcosh", method = "C")
  })
  print(time_base)
}

# Detailed timing breakdown
cat("\n=== Detailed timing breakdown ===\n")

# Time just the whitening step
cat("\nWhitening step:\n")
X_centered <- scale(t(X_expression), center = TRUE, scale = FALSE)
system.time({
  result_whiten <- RcppICA:::cpp_whiten(X_centered, n_comp, 0L)
})

# Time just the C++ call
cat("\nC++ backend only:\n")
system.time({
  result_cpp <- RcppICA:::cpp_fastICA_dense(
    X = X_centered,
    n_comp = as.integer(n_comp),
    alg_type = 0L,  # parallel
    fun_type = 0L,  # logcosh
    alpha = 1.0,
    whiten_method = 0L,  # SVD
    max_iter = 100L,
    tol = 1e-4,
    verbose = FALSE,
    n_threads = 1L,
    seed = 42L
  )
})

# Repeated runs for accuracy
cat("\n=== Repeated runs (5 times) ===\n")
times_rcpp <- numeric(5)
times_base <- numeric(5)

for (i in 1:5) {
  times_rcpp[i] <- system.time({
    fastICA(X_expression, n.comp = n_comp, maxit = 100,
            alg.typ = "parallel", n.threads = 1, seed = 42 + i)
  })[["elapsed"]]

  if (requireNamespace("fastICA", quietly = TRUE)) {
    times_base[i] <- system.time({
      fastICA::fastICA(X_expression, n.comp = n_comp, maxit = 100,
                       fun = "logcosh", method = "C")
    })[["elapsed"]]
  }
}

cat("RcppICA times:", round(times_rcpp, 3), "\n")
cat("RcppICA mean:", round(mean(times_rcpp), 3), "sd:", round(sd(times_rcpp), 3), "\n")
if (requireNamespace("fastICA", quietly = TRUE)) {
  cat("Base fastICA times:", round(times_base, 3), "\n")
  cat("Base fastICA mean:", round(mean(times_base), 3), "sd:", round(sd(times_base), 3), "\n")
  cat("Speedup ratio:", round(mean(times_base) / mean(times_rcpp), 2), "x\n")
}

# Check if results are similar
cat("\n=== Checking correctness ===\n")
if (exists("result_base")) {
  S_rcpp <- t(result_rcpp@S)
  S_base <- result_base$S

  cors <- abs(cor(S_rcpp, S_base))
  best_match <- apply(cors, 2, max)
  cat("Source recovery correlations:", round(best_match, 3), "\n")
  cat("Mean correlation:", round(mean(best_match), 3), "\n")
}

cat("\nRcppICA converged:", result_rcpp@misc$converged, "\n")
cat("RcppICA iterations:", result_rcpp@misc$iterations, "\n")
if (exists("result_base")) {
  cat("Base fastICA iterations:", result_base$iterations, "\n")
}
