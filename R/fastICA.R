#' Fast Independent Component Analysis
#'
#' Performs Independent Component Analysis using the FastICA algorithm
#' with optimized C++ implementation via Eigen and OpenMP.
#'
#' @param X A numeric matrix (observations x variables) or sparse dgCMatrix.
#'   Rows are observations, columns are variables.
#' @param n.comp Number of independent components to extract.
#'   Must be less than or equal to min(nrow(X), ncol(X)).
#' @param alg.typ Algorithm type: \code{"parallel"} (default, symmetric
#'   orthogonalization) or \code{"deflation"} (sequential extraction).
#' @param fun Nonlinearity function: \code{"logcosh"} (default, most robust),
#'   \code{"exp"} (good for super-Gaussian), or \code{"cube"} (kurtosis-based).
#' @param alpha Parameter for logcosh nonlinearity, between 1 and 2.
#'   Default is 1.0. Higher values approach kurtosis.
#' @param whiten.method Whitening method: \code{"spectra"} (default, fastest -
#'   computes only top-k eigenvalues using Lanczos iteration), \code{"eigen"}
#'   (full eigendecomposition of covariance matrix), or \code{"svd"}
#'   (SVD-based, most numerically stable for ill-conditioned data).
#' @param maxit Maximum number of iterations. Default is 200.
#' @param tol Convergence tolerance. Algorithm stops when the change in
#'   unmixing matrix is below this threshold. Default is 1e-4.
#' @param verbose Logical, print progress messages. Default is FALSE.
#' @param n.threads Number of OpenMP threads for parallel algorithm.
#'   Default is 0 (auto-detect maximum available).
#' @param seed Random seed for reproducibility. Default is NULL (random).
#'
#' @return An \code{\linkS4class{ICAResult}} object containing:
#' \describe{
#'   \item{S}{Matrix of independent components (n.comp x n observations).
#'     Each row is one independent component.}
#'   \item{A}{Mixing matrix (variables x n.comp). Maps from IC space to
#'     original variable space: X_centered = A \%*\% S}
#'   \item{W}{Unmixing matrix (n.comp x whitened dimension). Applied to
#'     whitened data to get ICs: S = W \%*\% X_whitened}
#'   \item{K}{Whitening matrix (variables x whitened dimension). Projects
#'     centered data to whitened space: X_whitened = X_centered \%*\% K}
#'   \item{center}{Column means used for centering}
#'   \item{iterations}{Number of iterations performed}
#'   \item{converged}{Logical indicating if algorithm converged}
#' }
#'
#' @details
#' FastICA is a fixed-point iteration algorithm that maximizes non-Gaussianity
#' to find independent components. This implementation uses:
#'
#' \itemize{
#'   \item Eigen C++ library for fast linear algebra
#'   \item OpenMP for parallel computation in the symmetric algorithm
#'   \item Spectra library for truncated eigendecomposition (default, fastest)
#'   \item Pre-allocated matrices to minimize memory allocation overhead
#' }
#'
#' The \strong{parallel} algorithm extracts all components simultaneously
#' using symmetric orthogonalization. It is typically faster and more stable.
#'
#' The \strong{deflation} algorithm extracts components one at a time using
#' Gram-Schmidt orthogonalization. It may be preferred when extracting only
#' a few components.
#'
#' @section Nonlinearity Functions:
#' \describe{
#'   \item{logcosh}{g(u) = tanh(alpha*u). Most robust choice, works well for
#'     both sub- and super-Gaussian distributions.}
#'   \item{exp}{g(u) = u*exp(-u^2/2). Good for super-Gaussian distributions
#'     with heavy tails.}
#'   \item{cube}{g(u) = u^3. Simplest, equivalent to maximizing kurtosis.
#'     Less robust to outliers.}
#' }
#'
#' @examples
#' # Generate mixed signals (cocktail party problem)
#' set.seed(42)
#' n <- 1000
#'
#' # Create two independent source signals
#' S <- cbind(
#'   sin(seq(0, 8*pi, length.out = n)),           # Sinusoid
#'   ((1:n %% 50) - 25) / 25                       # Sawtooth
#' )
#'
#' # Create random mixing matrix
#' A <- matrix(c(0.5, 0.8, 0.6, 0.3), 2, 2)
#'
#' # Mix the signals
#' X <- S %*% A
#'
#' # Recover independent components
#' result <- fastICA(X, n.comp = 2)
#'
#' # Check convergence
#' print(result)
#'
#' # Compare recovered vs original (allowing for sign/order permutation)
#' cor(t(result@@S), S)
#'
#' @seealso \code{\link{predict,ICAResult-method}} for projecting new data
#'
#' @references
#' Hyvarinen, A. and Oja, E. (2000). Independent Component Analysis:
#' Algorithms and Applications. Neural Networks, 13(4-5):411-430.
#'
#' @export
fastICA <- function(X,
                    n.comp,
                    alg.typ = c("parallel", "deflation"),
                    fun = c("logcosh", "exp", "cube"),
                    alpha = 1.0,
                    whiten.method = c("spectra", "eigen", "svd"),
                    maxit = 200L,
                    tol = 1e-4,
                    verbose = FALSE,
                    n.threads = 0L,
                    seed = NULL) {

    # Match arguments
    alg.typ <- match.arg(alg.typ)
    fun <- match.arg(fun)
    whiten.method <- match.arg(whiten.method)

    # Input validation
    is_sparse <- inherits(X, "dgCMatrix")

    if (!is.matrix(X) && !is_sparse) {
        X <- as.matrix(X)
    }

    if (!is_sparse && !is.numeric(X)) {
        stop("X must be a numeric matrix or dgCMatrix")
    }

    n <- nrow(X)
    m <- ncol(X)

    if (missing(n.comp)) {
        stop("n.comp must be specified")
    }

    n.comp <- as.integer(n.comp)
    if (n.comp < 1) {
        stop("n.comp must be at least 1")
    }
    if (n.comp > min(n, m)) {
        stop("n.comp cannot exceed min(nrow(X), ncol(X))")
    }

    # Validate alpha for logcosh
    if (fun == "logcosh" && (alpha < 1 || alpha > 2)) {
        stop("alpha must be between 1 and 2 for logcosh nonlinearity")
    }

    maxit <- as.integer(maxit)
    n.threads <- as.integer(n.threads)

    # Handle seed
    if (is.null(seed)) {
        seed <- 0L  # 0 means use random seed in C++
    } else {
        seed <- as.integer(seed)
    }

    # Map parameters to integer codes
    fun_code <- switch(fun, logcosh = 0L, exp = 1L, cube = 2L)
    alg_code <- switch(alg.typ, parallel = 0L, deflation = 1L)
    whiten_code <- switch(whiten.method, svd = 0L, eigen = 1L, spectra = 2L)

    # Route sparse and dense matrices to appropriate implementations
    if (is_sparse) {
        # Use sparse-aware implementation (major memory savings!)
        # This avoids materializing the dense centered matrix
        result <- cpp_fastICA_sparse(
            X_sparse = X,
            n_comp = n.comp,
            alg_type = alg_code,
            fun_type = fun_code,
            alpha = alpha,
            whiten_method = whiten_code,
            max_iter = maxit,
            tol = tol,
            verbose = verbose,
            n_threads = n.threads,
            seed = seed
        )
    } else {
        # Use standard dense implementation
        result <- cpp_fastICA_dense(
            X = X,
            n_comp = n.comp,
            alg_type = alg_code,
            fun_type = fun_code,
            alpha = alpha,
            whiten_method = whiten_code,
            max_iter = maxit,
            tol = tol,
            verbose = verbose,
            n_threads = n.threads,
            seed = seed
        )
    }

    # Construct S4 result object with misc metadata
    misc <- list(
        call = match.call(),
        iterations = as.integer(result$iterations),
        converged = as.logical(result$converged),
        algorithm = alg.typ,
        nonlinearity = fun,
        alpha = alpha,
        whiten_method = whiten.method,
        n_components = n.comp,
        tol = tol,
        maxit = maxit
    )

    new("ICAResult",
        S = result$S,
        A = result$A,
        W = result$W,
        K = result$K,
        center = as.numeric(result$center),
        misc = misc)
}
