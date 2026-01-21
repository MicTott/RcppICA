## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)
library(RcppICA)
set.seed(42)

## ----cocktail-party-----------------------------------------------------------
# Create two independent source signals
n <- 1000
time <- seq(0, 10, length.out = n)

# Source 1: Sinusoidal wave
source1 <- sin(2 * pi * time)

# Source 2: Sawtooth wave
source2 <- ((1:n %% 50) - 25) / 25

# Combine sources into matrix (observations × components)
S_true <- cbind(source1, source2)

# Create a random mixing matrix
A_true <- matrix(c(0.6, 0.8, 0.7, 0.4), 2, 2)

# Mix the signals
X <- S_true %*% A_true

# Plot mixed signals
par(mfrow = c(2, 1), mar = c(3, 4, 2, 1))
plot(time, X[, 1], type = "l", main = "Mixed Signal 1",
     ylab = "Amplitude", xlab = "")
plot(time, X[, 2], type = "l", main = "Mixed Signal 2",
     ylab = "Amplitude", xlab = "Time")

## ----separate-signals---------------------------------------------------------
# Run FastICA
result <- fastICA(X, n.comp = 2)

# Display results
print(result)

# Check how well we recovered the sources
# (Note: ICA can't recover exact scale or order, so we check correlation)
recovered_vs_true <- abs(cor(t(result@S), S_true))
print(recovered_vs_true)

## ----plot-recovered-----------------------------------------------------------
# Plot recovered components
par(mfrow = c(2, 1), mar = c(3, 4, 2, 1))
plot(time, result@S[1, ], type = "l", main = "Recovered Component 1",
     ylab = "Amplitude", xlab = "")
plot(time, result@S[2, ], type = "l", main = "Recovered Component 2",
     ylab = "Amplitude", xlab = "Time")

## ----result-structure---------------------------------------------------------
# Independent components (2 components × 1000 observations)
dim(result@S)

# Mixing matrix (2 variables × 2 components)
dim(result@A)

# Unmixing matrix (2 components × 2 whitened dimensions)
dim(result@W)

# Whitening matrix
dim(result@K)

# Convergence information
result@misc$iterations
result@misc$converged

## ----methods-demo-------------------------------------------------------------
# Summary statistics
summary(result)

# Extract components
S <- components(result)
head(S, 3)

# Extract mixing matrix
A <- mixing(result)
print(A)

## ----predict-demo-------------------------------------------------------------
# Create new mixed data
new_observations <- matrix(rnorm(20), 10, 2)
new_components <- predict(result, new_observations)

dim(new_components)

## ----image-example------------------------------------------------------------
# Create three simple image patterns
n_pixels <- 50
pattern1 <- outer(1:n_pixels, 1:n_pixels, function(x, y) sin(x/5) * cos(y/5))
pattern2 <- outer(1:n_pixels, 1:n_pixels, function(x, y)
  ifelse((x + y) %% 20 < 10, 1, -1))
pattern3 <- outer(1:n_pixels, 1:n_pixels, function(x, y) exp(-((x-25)^2 + (y-25)^2)/200))

# Flatten to vectors
S_images <- cbind(
  as.vector(pattern1),
  as.vector(pattern2),
  as.vector(pattern3)
)

# Create mixing matrix
A_mix <- matrix(c(0.5, 0.3, 0.7,
                  0.8, 0.4, 0.2,
                  0.3, 0.9, 0.4), 3, 3, byrow = TRUE)

# Mix images
X_mixed <- S_images %*% A_mix

# Run ICA to separate
result_images <- fastICA(X_mixed, n.comp = 3, maxit = 200)

cat("Converged:", result_images@misc$converged, "\n")
cat("Iterations:", result_images@misc$iterations, "\n")

# Check recovery (allowing for permutation and sign flip)
recovery_corr <- abs(cor(t(result_images@S), S_images))
cat("\nRecovery correlations:\n")
print(round(recovery_corr, 3))

## ----algorithm-types, eval=FALSE----------------------------------------------
# # Parallel (symmetric) - default, usually faster and more stable
# result_parallel <- fastICA(X, n.comp = 2, alg.typ = "parallel")
# 
# # Deflation (sequential) - extracts components one at a time
# result_deflation <- fastICA(X, n.comp = 2, alg.typ = "deflation")

## ----nonlinearities, eval=FALSE-----------------------------------------------
# # LogCosh (default) - robust, works well for most data
# result_logcosh <- fastICA(X, n.comp = 2, fun = "logcosh")
# 
# # Exp - good for super-Gaussian (heavy-tailed) distributions
# result_exp <- fastICA(X, n.comp = 2, fun = "exp")
# 
# # Cube - simple kurtosis-based, fast but less robust
# result_cube <- fastICA(X, n.comp = 2, fun = "cube")

## ----whitening, eval=FALSE----------------------------------------------------
# # SVD-based (default) - more numerically stable
# result_svd <- fastICA(X, n.comp = 2, whiten.method = "svd")
# 
# # Eigendecomposition - can be faster for n >> m
# result_eigen <- fastICA(X, n.comp = 2, whiten.method = "eigen")

## ----performance-tuning, eval=FALSE-------------------------------------------
# # Control convergence
# result <- fastICA(X, n.comp = 2,
#                   tol = 1e-6,      # Stricter tolerance
#                   maxit = 500)     # More iterations allowed
# 
# # Control parallelization
# result <- fastICA(X, n.comp = 2,
#                   n.threads = 4)   # Use 4 threads (0 = auto-detect)
# 
# # Reproducibility
# result <- fastICA(X, n.comp = 2,
#                   seed = 42)       # Set random seed

## ----benchmark, eval=FALSE----------------------------------------------------
# # Generate larger dataset
# X_large <- matrix(rnorm(10000), 1000, 10)
# 
# # Time RcppICA
# system.time({
#   result_rcpp <- fastICA(X_large, n.comp = 5, maxit = 200)
# })
# 
# # If you have the fastICA package installed:
# if (requireNamespace("fastICA", quietly = TRUE)) {
#   system.time({
#     result_base <- fastICA::fastICA(X_large, n.comp = 5, maxit = 200)
#   })
# }

## ----centering----------------------------------------------------------------
# Data is automatically centered
# The center is stored for predictions
result@center

## ----choosing-components, eval=FALSE------------------------------------------
# # Start with fewer components than variables
# ncol(X)  # Number of variables
# result <- fastICA(X, n.comp = min(5, ncol(X)))
# 
# # You can extract fewer components than the full rank

## ----check-convergence--------------------------------------------------------
if (!result@misc$converged) {
  warning("ICA did not converge. Try increasing maxit or adjusting tol.")
}

## ----troubleshooting, eval=FALSE----------------------------------------------
# # Try more iterations
# result <- fastICA(X, n.comp = 2, maxit = 500)
# 
# # Try looser tolerance
# result <- fastICA(X, n.comp = 2, tol = 1e-3)
# 
# # Try different nonlinearity
# result <- fastICA(X, n.comp = 2, fun = "exp")
# 
# # Try deflation algorithm
# result <- fastICA(X, n.comp = 2, alg.typ = "deflation")

## ----session-info-------------------------------------------------------------
sessionInfo()

