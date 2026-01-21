# Convergence behavior tests for RcppICA

test_that("ICA converges within reasonable iterations", {
  set.seed(11111)
  n <- 200
  S <- cbind(sin(seq(0, 8*pi, length.out = n)), rnorm(n))
  A <- matrix(c(0.6, 0.8, 0.7, 0.4), 2, 2)
  X <- S %*% A

  result <- fastICA(X, n.comp = 2, maxit = 200, tol = 1e-4)

  expect_true(result@misc$converged)
  expect_lt(result@misc$iterations, 200)
  expect_gt(result@misc$iterations, 0)
})

test_that("Stricter tolerance requires more iterations", {
  set.seed(22222)
  X <- matrix(rnorm(300), 100, 3)

  # Loose tolerance
  result_loose <- fastICA(X, n.comp = 2, maxit = 500, tol = 1e-3, seed = 123)

  # Strict tolerance
  result_strict <- fastICA(X, n.comp = 2, maxit = 500, tol = 1e-6, seed = 123)

  # Stricter tolerance should need more iterations (usually)
  # Note: Sometimes it converges quickly anyway, so we just check both converge
  expect_true(result_loose@misc$converged)
  expect_true(result_strict@misc$converged)
  expect_gte(result_strict@misc$iterations, result_loose@misc$iterations)
})

test_that("Max iterations limit is respected", {
  set.seed(33333)
  X <- matrix(rnorm(400), 100, 4)

  # Very low max iterations
  result <- fastICA(X, n.comp = 3, maxit = 5, tol = 1e-8)

  # Should stop at max iterations
  expect_lte(result@misc$iterations, 5)

  # May or may not converge with so few iterations
  if (!result@misc$converged) {
    expect_equal(result@misc$iterations, 5)
  }
})

test_that("Same seed produces identical results", {
  X <- matrix(rnorm(200), 100, 2)

  # Run multiple times with same seed
  result1 <- fastICA(X, n.comp = 2, seed = 42, maxit = 100)
  result2 <- fastICA(X, n.comp = 2, seed = 42, maxit = 100)
  result3 <- fastICA(X, n.comp = 2, seed = 42, maxit = 100)

  # All should be identical
  expect_equal(result1@S, result2@S)
  expect_equal(result2@S, result3@S)
  expect_equal(result1@misc$iterations, result2@misc$iterations)
  expect_equal(result2@misc$iterations, result3@misc$iterations)
})

test_that("Different seeds produce different results", {
  X <- matrix(rnorm(200), 100, 2)

  # Run with different seeds
  result1 <- fastICA(X, n.comp = 2, seed = 111, maxit = 100)
  result2 <- fastICA(X, n.comp = 2, seed = 222, maxit = 100)

  # Results should differ (different random initializations)
  expect_false(isTRUE(all.equal(result1@S, result2@S, tolerance = 1e-6)))

  # But both should converge to valid solutions
  expect_true(result1@misc$converged)
  expect_true(result2@misc$converged)

  # And both should recover the sources (allowing for permutation/sign)
  # Check that both have orthogonal unmixing matrices
  check_orthogonal <- function(W) {
    WWT <- W %*% t(W)
    all(abs(diag(WWT) - 1) < 0.01) && all(abs(WWT[upper.tri(WWT)]) < 0.01)
  }

  expect_true(check_orthogonal(result1@W))
  expect_true(check_orthogonal(result2@W))
})

test_that("Thread count parameter works correctly", {
  skip_on_cran()  # Skip on CRAN as thread behavior may vary

  set.seed(44444)
  X <- matrix(rnorm(500), 100, 5)

  # Test with different thread counts
  # Note: Results should be similar regardless of thread count

  result_1thread <- fastICA(X, n.comp = 3, n.threads = 1,
                             seed = 555, maxit = 100, alg.typ = "parallel")

  result_4threads <- fastICA(X, n.comp = 3, n.threads = 4,
                              seed = 555, maxit = 100, alg.typ = "parallel")

  result_auto <- fastICA(X, n.comp = 3, n.threads = 0,
                          seed = 555, maxit = 100, alg.typ = "parallel")

  # All should converge
  expect_true(result_1thread@misc$converged)
  expect_true(result_4threads@misc$converged)
  expect_true(result_auto@misc$converged)

  # Results should be very similar (same seed, same algorithm)
  # Allow small numerical differences due to parallel computation order
  cors <- abs(cor(t(result_1thread@S), t(result_4threads@S)))
  best_match <- apply(cors, 1, max)
  expect_true(all(best_match > 0.99))
})

test_that("Verbose mode doesn't crash", {
  set.seed(55555)
  X <- matrix(rnorm(200), 100, 2)

  # Just check that verbose mode runs without error
  # Output may or may not be captured depending on C++ implementation
  result <- fastICA(X, n.comp = 2, maxit = 50, verbose = TRUE)

  # Should still converge
  expect_s4_class(result, "ICAResult")
})

test_that("Convergence is stable across algorithm types", {
  set.seed(66666)
  n <- 200
  S <- cbind(
    sin(seq(0, 10*pi, length.out = n)),
    runif(n, -1, 1)
  )
  A <- matrix(c(0.6, 0.8, 0.7, 0.5), 2, 2)
  X <- S %*% A

  # Both algorithms should converge
  result_deflation <- fastICA(X, n.comp = 2, alg.typ = "deflation",
                               maxit = 200, tol = 1e-5)
  result_parallel <- fastICA(X, n.comp = 2, alg.typ = "parallel",
                              maxit = 200, tol = 1e-5)

  expect_true(result_deflation@misc$converged)
  expect_true(result_parallel@misc$converged)

  # Both should have reasonable iteration counts
  expect_lt(result_deflation@misc$iterations, 200)
  expect_lt(result_parallel@misc$iterations, 200)
})

test_that("Convergence with challenging data", {
  set.seed(77777)
  n <- 300

  # Create sources with different characteristics
  S <- cbind(
    sin(seq(0, 20*pi, length.out = n)),           # High frequency
    sign(sin(seq(0, 5*pi, length.out = n))),      # Square wave
    rnorm(n)                                       # Gaussian
  )

  A <- matrix(runif(9, 0.3, 0.7), 3, 3)
  X <- S %*% A

  result <- fastICA(X, n.comp = 3, maxit = 500, tol = 1e-5)

  # Should still converge, though may take more iterations
  expect_true(result@misc$converged)
  expect_lt(result@misc$iterations, 500)

  # Should still recover sources reasonably well
  cors <- abs(cor(t(result@S), S))
  best_matches <- apply(cors, 1, max)
  expect_true(mean(best_matches) > 0.75)
})

test_that("Early convergence is detected", {
  set.seed(88888)
  n <- 200

  # Create very clear independent sources
  S <- cbind(
    rep(c(-1, 1), each = n/2),              # Step function
    rep(c(-1, 0, 1, 0), each = n/4)         # Different pattern
  )

  A <- matrix(c(0.7, 0.3, 0.3, 0.7), 2, 2)
  X <- S %*% A

  result <- fastICA(X, n.comp = 2, maxit = 500, tol = 1e-4)

  # Should converge quickly with such clear sources
  expect_true(result@misc$converged)
  expect_lt(result@misc$iterations, 50)
})

test_that("Non-convergence is properly reported", {
  set.seed(99999)
  X <- matrix(rnorm(200), 100, 2)

  # Force non-convergence with very strict tolerance and few iterations
  result <- fastICA(X, n.comp = 2, maxit = 2, tol = 1e-12)

  # Should complete without error
  expect_s4_class(result, "ICAResult")

  # Should report non-convergence (or happened to converge in 2 iterations)
  if (!result@misc$converged) {
    expect_equal(result@misc$iterations, 2)
  }
})

test_that("Convergence metrics are stored correctly", {
  set.seed(10101)
  X <- matrix(rnorm(200), 100, 2)

  result <- fastICA(X, n.comp = 2, maxit = 100, tol = 1e-4)

  # Check misc list has all expected convergence info
  expect_true("iterations" %in% names(result@misc))
  expect_true("converged" %in% names(result@misc))
  expect_true("tol" %in% names(result@misc))
  expect_true("maxit" %in% names(result@misc))

  # Values should be sensible
  expect_true(is.integer(result@misc$iterations))
  expect_true(is.logical(result@misc$converged))
  expect_equal(result@misc$tol, 1e-4)
  expect_equal(result@misc$maxit, 100)
})
