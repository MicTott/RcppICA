## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)
library(RcppICA)
set.seed(123)

## ----parallel-example---------------------------------------------------------
# Generate test data
n <- 500
S <- cbind(
  sin(seq(0, 10*pi, length.out = n)),
  sign(sin(seq(0, 5*pi, length.out = n))),
  rnorm(n)
)
A <- matrix(runif(9, 0.3, 0.7), 3, 3)
X <- S %*% A

# Parallel algorithm (default)
result_parallel <- fastICA(X, n.comp = 3, alg.typ = "parallel", maxit = 200)
cat("Parallel - Iterations:", result_parallel@misc$iterations,
    "Converged:", result_parallel@misc$converged, "\n")

## ----deflation-example--------------------------------------------------------
# Deflation algorithm
result_deflation <- fastICA(X, n.comp = 3, alg.typ = "deflation", maxit = 200)
cat("Deflation - Iterations:", result_deflation@misc$iterations,
    "Converged:", result_deflation@misc$converged, "\n")

## ----compare-algorithms-------------------------------------------------------
# Compare recovered sources (allowing for permutation)
cors_parallel <- abs(cor(t(result_parallel@S), S))
cors_deflation <- abs(cor(t(result_deflation@S), S))

cat("Parallel recovery (best match per component):\n")
print(round(apply(cors_parallel, 1, max), 3))

cat("\nDeflation recovery (best match per component):\n")
print(round(apply(cors_deflation, 1, max), 3))

## ----logcosh------------------------------------------------------------------
# LogCosh with different alpha values
result_alpha1 <- fastICA(X, n.comp = 3, fun = "logcosh", alpha = 1.0, maxit = 200)
result_alpha2 <- fastICA(X, n.comp = 3, fun = "logcosh", alpha = 2.0, maxit = 200)

cat("Alpha = 1.0 iterations:", result_alpha1@misc$iterations, "\n")
cat("Alpha = 2.0 iterations:", result_alpha2@misc$iterations, "\n")

## ----exp-function-------------------------------------------------------------
# Create super-Gaussian source (sparse)
S_super <- cbind(
  rt(n, df = 3),  # Heavy-tailed
  rcauchy(n, scale = 0.5),
  rnorm(n)
)
X_super <- S_super %*% A

# Exponential nonlinearity
result_exp <- fastICA(X_super, n.comp = 3, fun = "exp", maxit = 200)
cat("Exp function - Iterations:", result_exp@misc$iterations,
    "Converged:", result_exp@misc$converged, "\n")

# Compare with logcosh
result_log <- fastICA(X_super, n.comp = 3, fun = "logcosh", maxit = 200)
cat("LogCosh - Iterations:", result_log@misc$iterations, "\n")

## ----cube-function------------------------------------------------------------
# Cube nonlinearity
result_cube <- fastICA(X, n.comp = 3, fun = "cube", maxit = 200)
cat("Cube function - Iterations:", result_cube@misc$iterations,
    "Converged:", result_cube@misc$converged, "\n")

## ----svd-whitening------------------------------------------------------------
# SVD whitening (default)
result_svd <- fastICA(X, n.comp = 3, whiten.method = "svd", maxit = 200)
cat("SVD whitening - Converged:", result_svd@misc$converged, "\n")

## ----eigen-whitening----------------------------------------------------------
# Eigen whitening
result_eigen <- fastICA(X, n.comp = 3, whiten.method = "eigen", maxit = 200)
cat("Eigen whitening - Converged:", result_eigen@misc$converged, "\n")

# Both should give similar results
cors <- abs(cor(t(result_svd@S), t(result_eigen@S)))
cat("\nCorrelation between SVD and Eigen results:\n")
print(round(apply(cors, 1, max), 3))

## ----tolerance----------------------------------------------------------------
# Loose tolerance (faster, less accurate)
result_loose <- fastICA(X, n.comp = 3, tol = 1e-3, maxit = 500)

# Strict tolerance (slower, more accurate)
result_strict <- fastICA(X, n.comp = 3, tol = 1e-7, maxit = 500)

cat("Loose (1e-3) - Iterations:", result_loose@misc$iterations, "\n")
cat("Strict (1e-7) - Iterations:", result_strict@misc$iterations, "\n")

## ----max-iterations-----------------------------------------------------------
# May not converge with few iterations
result_few <- fastICA(X, n.comp = 3, maxit = 10)
cat("maxit=10 - Converged:", result_few@misc$converged,
    "Iterations:", result_few@misc$iterations, "\n")

# More iterations = more chances to converge
result_many <- fastICA(X, n.comp = 3, maxit = 500)
cat("maxit=500 - Converged:", result_many@misc$converged,
    "Iterations:", result_many@misc$iterations, "\n")

## ----parallelization, eval=FALSE----------------------------------------------
# # Auto-detect cores (default)
# result_auto <- fastICA(X, n.comp = 3, alg.typ = "parallel", n.threads = 0)
# 
# # Explicit thread count
# result_4threads <- fastICA(X, n.comp = 3, alg.typ = "parallel", n.threads = 4)
# 
# # Single-threaded (for comparison)
# result_1thread <- fastICA(X, n.comp = 3, alg.typ = "parallel", n.threads = 1)

## ----reproducibility----------------------------------------------------------
# Same seed = identical results
result1 <- fastICA(X, n.comp = 3, seed = 42, maxit = 200)
result2 <- fastICA(X, n.comp = 3, seed = 42, maxit = 200)

# Check they're identical
all.equal(result1@S, result2@S)

# Different seeds = different (but equivalent) results
result3 <- fastICA(X, n.comp = 3, seed = 999, maxit = 200)

# Different initialization, but should still recover sources
cors_12 <- abs(cor(t(result1@S), t(result2@S)))
cors_13 <- abs(cor(t(result1@S), t(result3@S)))

cat("Same seed correlation:\n")
print(round(apply(cors_12, 1, max), 3))

cat("\nDifferent seed correlation:\n")
print(round(apply(cors_13, 1, max), 3))

## ----audio-separation, eval=FALSE---------------------------------------------
# result <- fastICA(mixed_audio,
#                   n.comp = n_sources,
#                   alg.typ = "parallel",    # Fast, multiple sources
#                   fun = "logcosh",         # Robust default
#                   whiten.method = "svd",   # Stable
#                   maxit = 200)

## ----financial, eval=FALSE----------------------------------------------------
# result <- fastICA(returns,
#                   n.comp = n_factors,
#                   alg.typ = "parallel",
#                   fun = "exp",             # Super-Gaussian
#                   whiten.method = "svd",
#                   tol = 1e-6,             # High precision
#                   maxit = 500)

## ----image-processing, eval=FALSE---------------------------------------------
# result <- fastICA(image_pixels,
#                   n.comp = n_features,
#                   alg.typ = "parallel",
#                   fun = "logcosh",         # Robust
#                   alpha = 1.0,
#                   whiten.method = "svd",
#                   maxit = 300)

## ----exploratory, eval=FALSE--------------------------------------------------
# result <- fastICA(data,
#                   n.comp = 5,
#                   alg.typ = "parallel",
#                   fun = "logcosh",
#                   tol = 1e-3,             # Fast convergence
#                   maxit = 100)            # Fewer iterations

## ----research, eval=FALSE-----------------------------------------------------
# result <- fastICA(data,
#                   n.comp = n,
#                   alg.typ = "parallel",
#                   fun = "logcosh",
#                   alpha = 1.0,
#                   whiten.method = "svd",
#                   tol = 1e-8,             # Very strict
#                   maxit = 1000,           # Many iterations
#                   seed = 12345)           # Reproducible

## ----session-info-------------------------------------------------------------
sessionInfo()

