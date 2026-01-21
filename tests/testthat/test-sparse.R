test_that("Sparse matrix ICA works correctly", {
    library(Matrix)

    # Generate test data
    set.seed(42)
    n <- 500
    m <- 100

    # Create sparse data (10% density, typical for scRNA-seq)
    nnz <- round(n * m * 0.10)
    i <- sample(1:n, nnz, replace = TRUE)
    j <- sample(1:m, nnz, replace = TRUE)
    x <- rexp(nnz, rate = 0.1)

    X_sparse <- sparseMatrix(i = i, j = j, x = x, dims = c(n, m))
    X_dense <- as.matrix(X_sparse)

    # Run ICA on both sparse and dense
    result_sparse <- fastICA(X_sparse, n.comp = 5, seed = 123)
    result_dense <- fastICA(X_dense, n.comp = 5, seed = 123)

    # Results should be very similar (allowing for numerical differences)
    # Check dimensions
    expect_equal(dim(result_sparse@S), c(5, n))
    expect_equal(dim(result_sparse@A), c(m, 5))

    # Check convergence
    expect_true(result_sparse@misc$converged)
    expect_true(result_dense@misc$converged)

    # Check that results are highly correlated
    # Note: ICA allows for sign/order permutation, so we check max correlation
    cors <- abs(cor(t(result_sparse@S), t(result_dense@S)))
    expect_true(all(apply(cors, 1, max) > 0.95))
})

test_that("Sparse whitening matches dense whitening", {
    library(Matrix)

    set.seed(42)
    n <- 200
    m <- 50

    # Create sparse data
    nnz <- round(n * m * 0.15)
    i <- sample(1:n, nnz, replace = TRUE)
    j <- sample(1:m, nnz, replace = TRUE)
    x <- rnorm(nnz)

    X_sparse <- sparseMatrix(i = i, j = j, x = x, dims = c(n, m))
    X_dense <- as.matrix(X_sparse)

    # Whiten using internal function
    result_sparse <- RcppICA:::cpp_whiten_sparse(X_sparse, n_comp = 10, method = 2L)
    result_dense <- RcppICA:::cpp_whiten(X_dense, n_comp = 10, method = 2L)

    # Whitening matrices should be very similar
    K_cor <- abs(cor(as.vector(result_sparse$K), as.vector(result_dense$K)))
    expect_true(K_cor > 0.99)

    # Means should match
    expect_true(cor(result_sparse$mean, result_dense$mean) > 0.99)

    # Whitened data should be similar
    Xw_cor <- cor(as.vector(result_sparse$X_whitened), as.vector(result_dense$X_whitened))
    expect_true(abs(Xw_cor) > 0.95)
})

test_that("Sparse implementation handles edge cases", {
    library(Matrix)

    # Very sparse matrix (1% density)
    set.seed(123)
    n <- 1000
    m <- 200
    nnz <- round(n * m * 0.01)

    i <- sample(1:n, nnz, replace = TRUE)
    j <- sample(1:m, nnz, replace = TRUE)
    x <- rexp(nnz)

    X_sparse <- sparseMatrix(i = i, j = j, x = x, dims = c(n, m))

    # Should run without errors
    expect_no_error({
        result <- fastICA(X_sparse, n.comp = 5, seed = 42)
    })

    result <- fastICA(X_sparse, n.comp = 5, seed = 42)
    expect_equal(dim(result@S), c(5, n))
    expect_true(result@misc$converged)
})

test_that("Sparse ICA with different algorithms", {
    library(Matrix)

    set.seed(42)
    n <- 300
    m <- 80
    nnz <- round(n * m * 0.10)

    i <- sample(1:n, nnz, replace = TRUE)
    j <- sample(1:m, nnz, replace = TRUE)
    x <- rexp(nnz)

    X_sparse <- sparseMatrix(i = i, j = j, x = x, dims = c(n, m))

    # Parallel algorithm
    result_parallel <- fastICA(X_sparse, n.comp = 3, alg.typ = "parallel", seed = 42)
    expect_true(result_parallel@misc$converged)
    expect_equal(dim(result_parallel@S), c(3, n))

    # Deflation algorithm
    result_deflation <- fastICA(X_sparse, n.comp = 3, alg.typ = "deflation", seed = 42)
    expect_true(result_deflation@misc$converged)
    expect_equal(dim(result_deflation@S), c(3, n))

    # Results should be similar (allowing for algorithm differences)
    cors <- abs(cor(t(result_parallel@S), t(result_deflation@S)))
    expect_true(all(apply(cors, 1, max) > 0.85))
})

test_that("Sparse ICA with different nonlinearities", {
    library(Matrix)

    set.seed(42)
    n <- 250
    m <- 60
    nnz <- round(n * m * 0.12)

    i <- sample(1:n, nnz, replace = TRUE)
    j <- sample(1:m, nnz, replace = TRUE)
    x <- rnorm(nnz)

    X_sparse <- sparseMatrix(i = i, j = j, x = x, dims = c(n, m))

    # Test different nonlinearities
    result_logcosh <- fastICA(X_sparse, n.comp = 4, fun = "logcosh", seed = 42)
    result_exp <- fastICA(X_sparse, n.comp = 4, fun = "exp", seed = 42)
    result_cube <- fastICA(X_sparse, n.comp = 4, fun = "cube", seed = 42)

    # All should converge
    expect_true(result_logcosh@misc$converged)
    expect_true(result_exp@misc$converged)
    expect_true(result_cube@misc$converged)

    # All should have correct dimensions
    expect_equal(dim(result_logcosh@S), c(4, n))
    expect_equal(dim(result_exp@S), c(4, n))
    expect_equal(dim(result_cube@S), c(4, n))
})
