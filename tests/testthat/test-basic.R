# Basic functionality tests for RcppICA

test_that("fastICA runs on simple data", {
  set.seed(123)
  n <- 100

  # Create two independent sources
  S <- cbind(
    sin(seq(0, 8*pi, length.out = n)),
    ((1:n %% 20) - 10) / 10
  )

  # Create mixing matrix
  A <- matrix(c(0.5, 0.8, 0.6, 0.3), 2, 2)

  # Mix signals
  X <- S %*% A

  # Run ICA
  result <- fastICA(X, n.comp = 2, maxit = 100)

  # Basic checks
  expect_s4_class(result, "ICAResult")
  expect_equal(nrow(result@S), 2)
  expect_equal(ncol(result@S), n)
  expect_equal(nrow(result@A), 2)
  expect_equal(ncol(result@A), 2)
  expect_equal(nrow(result@W), 2)
  expect_equal(ncol(result@W), 2)
  expect_equal(nrow(result@K), 2)
  expect_equal(ncol(result@K), 2)
  expect_length(result@center, 2)
  expect_true(is.list(result@misc))
  expect_true(result@misc$converged)
})

test_that("ICAResult object has correct structure", {
  set.seed(456)
  X <- matrix(rnorm(200), 100, 2)
  result <- fastICA(X, n.comp = 2, maxit = 50)

  # Check slots exist
  expect_true(all(c("S", "A", "W", "K", "center", "misc") %in% slotNames(result)))

  # Check misc metadata
  expect_true("iterations" %in% names(result@misc))
  expect_true("converged" %in% names(result@misc))
  expect_true("call" %in% names(result@misc))
  expect_true("algorithm" %in% names(result@misc))
  expect_true("nonlinearity" %in% names(result@misc))
})

test_that("fastICA handles different dimensions", {
  set.seed(789)

  # n > p (more observations than variables)
  X1 <- matrix(rnorm(500), 100, 5)
  result1 <- fastICA(X1, n.comp = 3, maxit = 50)
  expect_equal(nrow(result1@S), 3)
  expect_equal(ncol(result1@S), 100)

  # n < p (fewer observations than variables) - should still work with reduced components
  X2 <- matrix(rnorm(200), 20, 10)
  result2 <- fastICA(X2, n.comp = 5, maxit = 50)
  expect_equal(nrow(result2@S), 5)
  expect_equal(ncol(result2@S), 20)

  # Square matrix
  X3 <- matrix(rnorm(400), 20, 20)
  result3 <- fastICA(X3, n.comp = 10, maxit = 50)
  expect_equal(nrow(result3@S), 10)
})

test_that("fastICA validates input parameters", {
  X <- matrix(rnorm(200), 100, 2)

  # Invalid n.comp
  expect_error(fastICA(X, n.comp = 0), "must be at least 1")
  expect_error(fastICA(X, n.comp = -1), "must be at least 1")

  # Invalid algorithm type
  expect_error(fastICA(X, n.comp = 2, alg.typ = "invalid"), "should be one of")

  # Invalid nonlinearity
  expect_error(fastICA(X, n.comp = 2, fun = "invalid"), "should be one of")

  # Invalid alpha
  expect_error(fastICA(X, n.comp = 2, alpha = 0.5), "between 1 and 2")
  expect_error(fastICA(X, n.comp = 2, alpha = 3), "between 1 and 2")

  # Invalid whitening method
  expect_error(fastICA(X, n.comp = 2, whiten.method = "invalid"), "should be one of")
})

test_that("fastICA handles edge cases", {
  set.seed(101)

  # Single component
  X <- matrix(rnorm(200), 100, 2)
  result <- fastICA(X, n.comp = 1, maxit = 50)
  expect_equal(nrow(result@S), 1)
  expect_equal(ncol(result@S), 100)

  # Maximum components (extract all)
  X2 <- matrix(rnorm(120), 40, 3)
  result2 <- fastICA(X2, n.comp = 3, maxit = 50)
  expect_equal(nrow(result2@S), 3)
})

test_that("show method works", {
  set.seed(202)
  X <- matrix(rnorm(200), 100, 2)
  result <- fastICA(X, n.comp = 2, maxit = 50)

  # Capture output
  output <- capture.output(show(result))

  # Check key information is displayed
  expect_true(any(grepl("ICA Result", output)))
  expect_true(any(grepl("Components:", output)))
  expect_true(any(grepl("Observations:", output)))
  expect_true(any(grepl("Variables:", output)))
  expect_true(any(grepl("Iterations:", output)))
  expect_true(any(grepl("Converged:", output)))
})

test_that("summary method works", {
  set.seed(303)
  X <- matrix(rnorm(200), 100, 2)
  result <- fastICA(X, n.comp = 2, maxit = 50)

  # Capture output
  output <- capture.output(summary(result))

  # Check detailed information is displayed
  expect_true(any(grepl("Independent Component Analysis", output)))
  expect_true(any(grepl("Dimensions:", output)))
  expect_true(any(grepl("Convergence:", output)))
  expect_true(any(grepl("Component Statistics:", output)))
})

test_that("predict method works", {
  set.seed(404)
  X_train <- matrix(rnorm(200), 100, 2)
  X_test <- matrix(rnorm(40), 20, 2)

  result <- fastICA(X_train, n.comp = 2, maxit = 50)

  # Predict on new data
  S_new <- predict(result, X_test)

  expect_true(is.matrix(S_new))
  expect_equal(nrow(S_new), 20)
  expect_equal(ncol(S_new), 2)

  # Predict without new data (return transposed S)
  S_orig <- predict(result)
  expect_equal(nrow(S_orig), 100)
  expect_equal(ncol(S_orig), 2)
})

test_that("components and mixing extractors work", {
  set.seed(505)
  X <- matrix(rnorm(200), 100, 2)
  result <- fastICA(X, n.comp = 2, maxit = 50)

  # Extract components
  S <- components(result)
  expect_identical(S, result@S)

  # Extract mixing matrix
  A <- mixing(result)
  expect_identical(A, result@A)
})

test_that("plot method works without errors", {
  set.seed(606)
  X <- matrix(rnorm(200), 100, 2)
  result <- fastICA(X, n.comp = 2, maxit = 50)

  # Plot components
  expect_silent({
    pdf(NULL)  # Don't actually display
    plot(result, type = "components")
    dev.off()
  })

  # Plot mixing matrix
  expect_silent({
    pdf(NULL)
    plot(result, type = "mixing")
    dev.off()
  })

  # Plot specific components
  expect_silent({
    pdf(NULL)
    plot(result, components = 1)
    dev.off()
  })
})

test_that("fastICA is reproducible with seed", {
  X <- matrix(rnorm(200), 100, 2)

  # Run with same seed
  result1 <- fastICA(X, n.comp = 2, seed = 12345, maxit = 50)
  result2 <- fastICA(X, n.comp = 2, seed = 12345, maxit = 50)

  # Results should be identical
  expect_equal(result1@S, result2@S)
  expect_equal(result1@A, result2@A)
  expect_equal(result1@W, result2@W)
  expect_equal(result1@misc$iterations, result2@misc$iterations)
})

test_that("fastICA handles constant columns", {
  set.seed(707)
  X <- matrix(rnorm(200), 100, 2)
  X[, 2] <- 1  # Constant column

  # May produce degenerate results but should not crash
  # The algorithm will likely not converge properly
  result <- fastICA(X, n.comp = 1, maxit = 50)
  expect_s4_class(result, "ICAResult")
})
