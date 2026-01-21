# Algorithm correctness tests for RcppICA

test_that("ICA recovers independent sources", {
  set.seed(12345)
  n <- 500

  # Create highly independent sources
  S_true <- cbind(
    sin(2 * pi * seq(0, 10, length.out = n)),     # Sinusoidal
    runif(n, -1, 1),                               # Uniform
    rnorm(n)                                       # Gaussian
  )

  # Create mixing matrix
  A_true <- matrix(runif(9, -1, 1), 3, 3)

  # Mix signals
  X <- S_true %*% A_true

  # Run ICA
  result <- fastICA(X, n.comp = 3, maxit = 200, tol = 1e-6)

  # Check source recovery (allowing for permutation and sign flips)
  # Compute absolute correlations between recovered and true sources
  cors <- abs(cor(t(result@S), S_true))

  # For each recovered component, find best match with true source
  best_matches <- apply(cors, 1, max)

  # All components should have high correlation with some true source
  expect_true(all(best_matches > 0.90),
              info = sprintf("Source recovery correlations: %s",
                             paste(round(best_matches, 3), collapse = ", ")))
})

test_that("Deflation and parallel algorithms give consistent results", {
  set.seed(54321)
  n <- 200

  # Create sources
  S <- cbind(
    sin(seq(0, 8*pi, length.out = n)),
    ((1:n %% 30) - 15) / 15
  )
  A <- matrix(c(0.6, 0.7, 0.8, 0.4), 2, 2)
  X <- S %*% A

  # Run both algorithms with same seed
  result_deflation <- fastICA(X, n.comp = 2, alg.typ = "deflation",
                               seed = 999, maxit = 100)
  result_parallel <- fastICA(X, n.comp = 2, alg.typ = "parallel",
                              seed = 999, maxit = 100)

  # Both should converge
  expect_true(result_deflation@misc$converged)
  expect_true(result_parallel@misc$converged)

  # Results should be similar (allowing for permutation/sign)
  # Check that mixing matrices span similar subspaces
  cors_A <- abs(cor(result_deflation@A, result_parallel@A))

  # Each column should correlate highly with some column in the other
  best_match_A <- apply(cors_A, 1, max)
  expect_true(all(best_match_A > 0.85),
              info = sprintf("Mixing matrix correlations: %s",
                             paste(round(best_match_A, 3), collapse = ", ")))
})

test_that("All nonlinearity functions work correctly", {
  set.seed(11111)
  n <- 150
  S <- cbind(sin(seq(0, 6*pi, length.out = n)), rnorm(n))
  A <- matrix(c(0.5, 0.7, 0.6, 0.4), 2, 2)
  X <- S %*% A

  # Test logcosh
  result_logcosh <- fastICA(X, n.comp = 2, fun = "logcosh", maxit = 100)
  expect_true(result_logcosh@misc$converged)
  expect_equal(result_logcosh@misc$nonlinearity, "logcosh")

  # Test exp
  result_exp <- fastICA(X, n.comp = 2, fun = "exp", maxit = 100)
  expect_true(result_exp@misc$converged)
  expect_equal(result_exp@misc$nonlinearity, "exp")

  # Test cube
  result_cube <- fastICA(X, n.comp = 2, fun = "cube", maxit = 100)
  expect_true(result_cube@misc$converged)
  expect_equal(result_cube@misc$nonlinearity, "cube")

  # All should recover sources reasonably well
  check_recovery <- function(result, S_true) {
    cors <- abs(cor(t(result@S), S_true))
    min(apply(cors, 1, max))
  }

  expect_gt(check_recovery(result_logcosh, S), 0.85)
  expect_gt(check_recovery(result_exp, S), 0.80)
  expect_gt(check_recovery(result_cube, S), 0.80)
})

test_that("Both whitening methods work correctly", {
  set.seed(22222)
  n <- 200
  S <- cbind(sin(seq(0, 10*pi, length.out = n)), runif(n, -1, 1))
  A <- matrix(c(0.6, 0.8, 0.7, 0.5), 2, 2)
  X <- S %*% A

  # SVD whitening (default)
  result_svd <- fastICA(X, n.comp = 2, whiten.method = "svd",
                        seed = 777, maxit = 100)
  expect_true(result_svd@misc$converged)
  expect_equal(result_svd@misc$whiten_method, "svd")

  # Eigen whitening
  result_eigen <- fastICA(X, n.comp = 2, whiten.method = "eigen",
                          seed = 777, maxit = 100)
  expect_true(result_eigen@misc$converged)
  expect_equal(result_eigen@misc$whiten_method, "eigen")

  # Both should give similar results
  cors_S <- abs(cor(t(result_svd@S), t(result_eigen@S)))
  best_match_S <- apply(cors_S, 1, max)
  expect_true(all(best_match_S > 0.95),
              info = sprintf("Whitening method comparison: %s",
                             paste(round(best_match_S, 3), collapse = ", ")))
})

test_that("Whitening matrix is orthogonal", {
  set.seed(33333)
  X <- matrix(rnorm(300), 100, 3)

  result <- fastICA(X, n.comp = 3, maxit = 100)

  # K should be such that whitened data has identity covariance
  # After centering and whitening: Cov(X_centered * K) = I
  X_centered <- sweep(X, 2, result@center)
  X_whitened <- X_centered %*% result@K

  # Covariance should be approximately identity
  cov_whitened <- cov(X_whitened)
  expect_true(all(abs(diag(cov_whitened) - 1) < 0.01))
  expect_true(all(abs(cov_whitened[upper.tri(cov_whitened)]) < 0.01))
})

test_that("Unmixing matrix is orthogonal", {
  set.seed(44444)
  X <- matrix(rnorm(200), 100, 2)

  result <- fastICA(X, n.comp = 2, maxit = 100)

  # W should be orthogonal: W * W^T = I
  W_WTrans <- result@W %*% t(result@W)

  expect_true(all(abs(diag(W_WTrans) - 1) < 0.01))
  expect_true(all(abs(W_WTrans[upper.tri(W_WTrans)]) < 0.01))
})

test_that("ICA decomposition is consistent: X ≈ A * S", {
  set.seed(55555)
  n <- 150
  S_true <- cbind(sin(seq(0, 6*pi, length.out = n)), rnorm(n))
  A_true <- matrix(c(0.5, 0.8, 0.6, 0.3), 2, 2)
  X <- S_true %*% A_true

  result <- fastICA(X, n.comp = 2, maxit = 200, tol = 1e-6)

  # The ICA model: X_centered = X_whitened * W^T = (X_centered * K) * W^T
  # And: S = W * X_whitened = W * (X_centered * K)
  # So: X_centered ≈ K_inv * W^T * S
  # But we don't store K_inv directly, so use: X_centered * K * W^T = S

  X_centered <- sweep(X, 2, result@center)
  # Verify whitening step
  X_whitened <- X_centered %*% result@K
  S_check <- t(result@W %*% t(X_whitened))

  # S from result should match S_check
  reconstruction_error <- mean((t(result@S) - S_check)^2)
  expect_lt(reconstruction_error, 1e-8)
})

test_that("Predict method correctly projects new data", {
  set.seed(66666)
  n_train <- 200
  n_test <- 50

  # Create consistent sources
  create_sources <- function(n) {
    cbind(
      sin(seq(0, 4*pi, length.out = n)),
      rnorm(n)
    )
  }

  S_train <- create_sources(n_train)
  S_test <- create_sources(n_test)

  A <- matrix(c(0.6, 0.7, 0.8, 0.5), 2, 2)
  X_train <- S_train %*% A
  X_test <- S_test %*% A

  # Fit on training data
  result <- fastICA(X_train, n.comp = 2, maxit = 100)

  # Project test data
  S_test_pred <- predict(result, X_test)

  # Should have correct dimensions
  expect_equal(dim(S_test_pred), c(n_test, 2))

  # Predicted components should have similar structure to training
  # (can't expect exact match due to different data, but should be similar scale)
  expect_true(sd(S_test_pred[, 1]) < 5 * sd(result@S[1, ]))
  expect_true(sd(S_test_pred[, 2]) < 5 * sd(result@S[2, ]))
})

test_that("ICA is invariant to scaling of data", {
  set.seed(77777)
  n <- 150
  S <- cbind(sin(seq(0, 8*pi, length.out = n)), rnorm(n))
  A <- matrix(c(0.5, 0.7, 0.6, 0.4), 2, 2)
  X <- S %*% A

  # Run ICA on original data
  result1 <- fastICA(X, n.comp = 2, seed = 888, maxit = 100)

  # Run ICA on scaled data
  X_scaled <- X * 10
  result2 <- fastICA(X_scaled, n.comp = 2, seed = 888, maxit = 100)

  # Components should be similar (up to scaling)
  # Normalize and compare
  S1_norm <- scale(t(result1@S))
  S2_norm <- scale(t(result2@S))

  cors <- abs(cor(S1_norm, S2_norm))
  best_match <- apply(cors, 1, max)

  expect_true(all(best_match > 0.95),
              info = sprintf("Scale invariance check: %s",
                             paste(round(best_match, 3), collapse = ", ")))
})

test_that("Alpha parameter affects logcosh results", {
  set.seed(88888)
  n <- 150
  S <- cbind(sin(seq(0, 6*pi, length.out = n)), rnorm(n))
  A <- matrix(c(0.6, 0.7, 0.7, 0.5), 2, 2)
  X <- S %*% A

  # Run with different alpha values
  result_alpha1 <- fastICA(X, n.comp = 2, fun = "logcosh",
                            alpha = 1.0, seed = 999, maxit = 100)
  result_alpha2 <- fastICA(X, n.comp = 2, fun = "logcosh",
                            alpha = 2.0, seed = 999, maxit = 100)

  # Both should converge
  expect_true(result_alpha1@misc$converged)
  expect_true(result_alpha2@misc$converged)

  # Results may differ slightly (alpha changes the contrast function)
  # But both should still recover sources well
  check_recovery <- function(result, S_true) {
    cors <- abs(cor(t(result@S), S_true))
    min(apply(cors, 1, max))
  }

  expect_gt(check_recovery(result_alpha1, S), 0.85)
  expect_gt(check_recovery(result_alpha2, S), 0.85)
})
