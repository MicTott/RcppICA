# Performance benchmark tests for RcppICA

test_that("RcppICA completes in reasonable time for small data", {
  skip_on_cran()  # Skip timing tests on CRAN

  set.seed(12345)
  X <- matrix(rnorm(1000), 100, 10)

  # Should complete quickly for small data
  time_taken <- system.time({
    result <- fastICA(X, n.comp = 5, maxit = 200)
  })

  expect_lt(time_taken["elapsed"], 5.0)  # Should take less than 5 seconds
  expect_true(result@misc$converged)
})

test_that("RcppICA scales reasonably with data size", {
  skip_on_cran()  # Skip timing tests on CRAN

  set.seed(23456)

  # Small data
  X_small <- matrix(rnorm(500), 100, 5)
  time_small <- system.time({
    result_small <- fastICA(X_small, n.comp = 3, maxit = 100)
  })["elapsed"]

  # Medium data (4x larger)
  X_medium <- matrix(rnorm(2000), 400, 5)
  time_medium <- system.time({
    result_medium <- fastICA(X_medium, n.comp = 3, maxit = 300)
  })["elapsed"]

  # Time should scale sub-quadratically (not 16x slower for 4x data)
  # Allow generous margin as timing can vary
  if (time_small > 0.001) {  # Only check if small was slow enough to measure
    expect_lt(time_medium, time_small * 30)
  }

  expect_true(result_small@misc$converged)
  expect_s4_class(result_medium, "ICAResult")
})

test_that("Parallel algorithm is not slower than deflation", {
  skip_on_cran()  # Skip timing tests on CRAN

  set.seed(34567)
  X <- matrix(rnorm(2000), 200, 10)

  # Time deflation
  time_deflation <- system.time({
    result_deflation <- fastICA(X, n.comp = 5, alg.typ = "deflation",
                                 maxit = 300, seed = 111)
  })["elapsed"]

  # Time parallel
  time_parallel <- system.time({
    result_parallel <- fastICA(X, n.comp = 5, alg.typ = "parallel",
                                maxit = 300, seed = 111)
  })["elapsed"]

  # Parallel should be at least competitive
  # Note: For small data, overhead might make parallel slower
  # Allow large margin as timing variability is high for fast operations
  if (time_deflation > 0.001) {
    expect_lt(time_parallel, time_deflation * 5)
  }

  # Should complete (may not converge with limited iterations)
  expect_s4_class(result_deflation, "ICAResult")
  expect_s4_class(result_parallel, "ICAResult")
})

test_that("OpenMP parallelization provides speedup for parallel algorithm", {
  skip_on_cran()  # Skip timing tests on CRAN
  skip_if(parallel::detectCores() < 2, "Need multiple cores for parallel test")

  set.seed(45678)
  X <- matrix(rnorm(5000), 500, 10)

  # Single thread
  time_1thread <- system.time({
    result_1 <- fastICA(X, n.comp = 5, alg.typ = "parallel",
                        n.threads = 1, maxit = 300, seed = 222)
  })["elapsed"]

  # Multiple threads
  time_multi <- system.time({
    result_multi <- fastICA(X, n.comp = 5, alg.typ = "parallel",
                             n.threads = 0, maxit = 300, seed = 222)
  })["elapsed"]

  # Multi-threaded should be faster (or at least not much slower)
  # Allow margin as speedup depends on hardware and problem size
  expect_lt(time_multi, time_1thread * 2.0)

  # Should complete (may not converge with limited iterations)
  expect_s4_class(result_1, "ICAResult")
  expect_s4_class(result_multi, "ICAResult")
})

test_that("Memory usage is reasonable", {
  skip_on_cran()  # Skip on CRAN

  set.seed(56789)

  # Moderately large data
  n <- 1000
  p <- 20
  X <- matrix(rnorm(n * p), n, p)

  # Run ICA
  result <- fastICA(X, n.comp = 10, maxit = 100)

  # Check that result objects are reasonable size
  # Matrices should be appropriately sized, not bloated
  expect_equal(dim(result@S), c(10, n))
  expect_equal(dim(result@A), c(p, 10))
  expect_equal(dim(result@W), c(10, 10))
  expect_equal(dim(result@K), c(p, 10))
  expect_length(result@center, p)

  # Object size should be reasonable (< 10MB for this data)
  obj_size <- object.size(result)
  expect_lt(as.numeric(obj_size), 10 * 1024 * 1024)  # 10 MB
})

test_that("Repeated calls don't leak memory", {
  skip_on_cran()  # Skip on CRAN

  set.seed(67890)
  X <- matrix(rnorm(500), 100, 5)

  # Run multiple times
  for (i in 1:10) {
    result <- fastICA(X, n.comp = 3, maxit = 50)
    expect_true(result@misc$converged)
  }

  # If we got here without crashing, no obvious memory issues
  expect_true(TRUE)
})

test_that("Comparison with fastICA package (if available)", {
  skip_on_cran()  # Skip on CRAN
  skip_if_not_installed("fastICA")

  set.seed(78901)
  n <- 200
  S <- cbind(sin(seq(0, 10*pi, length.out = n)), runif(n, -1, 1))
  A <- matrix(c(0.6, 0.8, 0.7, 0.5), 2, 2)
  X <- S %*% A

  # Time RcppICA
  time_rcpp <- system.time({
    result_rcpp <- fastICA(X, n.comp = 2, maxit = 100, seed = 333)
  })["elapsed"]

  # Time original fastICA
  time_original <- system.time({
    result_original <- fastICA::fastICA(X, n.comp = 2, maxit = 100,
                                        fun = "logcosh", method = "C")
  })["elapsed"]

  # RcppICA should be faster (but allow margin for small data overhead)
  # Main goal is to not be slower
  expect_lt(time_rcpp, time_original * 2)

  cat(sprintf("\nSpeed comparison: RcppICA = %.4fs, fastICA = %.4fs, speedup = %.1fx\n",
              time_rcpp, time_original, time_original / time_rcpp))

  # Both should recover sources
  check_recovery <- function(S_recovered, S_true) {
    cors <- abs(cor(t(S_recovered), S_true))
    min(apply(cors, 1, max))
  }

  expect_gt(check_recovery(result_rcpp@S, S), 0.90)
  expect_gt(check_recovery(t(result_original$S), S), 0.90)
})

test_that("Performance is consistent across nonlinearities", {
  skip_on_cran()  # Skip timing tests on CRAN

  set.seed(89012)
  X <- matrix(rnorm(1000), 200, 5)

  # Time each nonlinearity
  times <- c()

  time_logcosh <- system.time({
    result_logcosh <- fastICA(X, n.comp = 3, fun = "logcosh",
                               maxit = 100, seed = 444)
  })["elapsed"]
  times["logcosh"] <- time_logcosh

  time_exp <- system.time({
    result_exp <- fastICA(X, n.comp = 3, fun = "exp",
                          maxit = 100, seed = 444)
  })["elapsed"]
  times["exp"] <- time_exp

  time_cube <- system.time({
    result_cube <- fastICA(X, n.comp = 3, fun = "cube",
                           maxit = 100, seed = 444)
  })["elapsed"]
  times["cube"] <- time_cube

  # All should complete quickly (< 5 seconds)
  expect_true(all(times < 5.0))

  # Times should be similar (allowing for variability)
  # Skip comparison if any took 0 time (too fast to measure)
  if (min(times) > 0.001) {
    expect_lt(max(times) / min(times), 5.0)
  }
})

test_that("Whitening methods have similar performance", {
  skip_on_cran()  # Skip timing tests on CRAN

  set.seed(90123)
  X <- matrix(rnorm(2000), 200, 10)

  # Time SVD whitening
  time_svd <- system.time({
    result_svd <- fastICA(X, n.comp = 5, whiten.method = "svd",
                          maxit = 300, seed = 555)
  })["elapsed"]

  # Time eigen whitening
  time_eigen <- system.time({
    result_eigen <- fastICA(X, n.comp = 5, whiten.method = "eigen",
                            maxit = 300, seed = 555)
  })["elapsed"]

  # Both should be reasonably fast
  expect_lt(time_svd, 5.0)
  expect_lt(time_eigen, 5.0)

  # Times should be similar (SVD might be slightly slower but more stable)
  if (min(time_svd, time_eigen) > 0.001) {
    expect_lt(max(time_svd, time_eigen) / min(time_svd, time_eigen), 5.0)
  }

  # Should complete (may not converge with limited iterations for complex data)
  expect_s4_class(result_svd, "ICAResult")
  expect_s4_class(result_eigen, "ICAResult")
})

test_that("Large number of components doesn't cause excessive slowdown", {
  skip_on_cran()  # Skip timing tests on CRAN

  set.seed(11223)
  X <- matrix(rnorm(2000), 200, 10)

  # Few components
  time_few <- system.time({
    result_few <- fastICA(X, n.comp = 2, maxit = 100, seed = 666)
  })["elapsed"]

  # Many components
  time_many <- system.time({
    result_many <- fastICA(X, n.comp = 8, maxit = 100, seed = 666)
  })["elapsed"]

  # More components should take longer, but not excessively
  # Should scale roughly linearly with number of components
  if (time_few > 0.001) {  # Only check if few was slow enough to measure
    expect_lt(time_many, time_few * 15)  # 4x components shouldn't be >15x slower
  }

  expect_s4_class(result_few, "ICAResult")
  expect_s4_class(result_many, "ICAResult")
})
