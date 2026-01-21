// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

#include "../inst/include/RcppICA.h"
#include <RcppEigen.h>

using namespace Rcpp;
using namespace RcppICA;

// Helper function to run ICA with a specific nonlinearity
template<typename Nonlinearity>
List runICA(const Eigen::MatrixXd& X,
            int n_comp,
            int alg_type,
            int whiten_method,
            int max_iter,
            double tol,
            bool verbose,
            int n_threads,
            unsigned int seed,
            Nonlinearity nonlin) {

    ICAResult<double> result;

    WhiteningMethod wmethod = static_cast<WhiteningMethod>(whiten_method);

    if (alg_type == 0) {
        // Parallel (symmetric) algorithm
        ICAParallel<double, Nonlinearity> ica(n_comp, max_iter, tol, verbose,
                                               n_threads, nonlin, wmethod);
        result = ica.fit(X, seed);
    } else {
        // Deflation algorithm
        ICADeflation<double, Nonlinearity> ica(n_comp, max_iter, tol, verbose,
                                                nonlin, wmethod);
        result = ica.fit(X, seed);
    }

    return List::create(
        Named("S") = result.S,
        Named("A") = result.A,
        Named("W") = result.W,
        Named("K") = result.K,
        Named("center") = result.center,
        Named("iterations") = result.iterations,
        Named("converged") = result.converged
    );
}

//' Fast Independent Component Analysis (C++ Implementation)
//'
//' @param X Numeric matrix (observations x variables)
//' @param n_comp Number of independent components to extract
//' @param alg_type Algorithm type: 0 = parallel (symmetric), 1 = deflation
//' @param fun_type Nonlinearity: 0 = logcosh, 1 = exp, 2 = cube
//' @param alpha Parameter for logcosh (between 1 and 2)
//' @param whiten_method Whitening: 0 = SVD, 1 = eigendecomposition
//' @param max_iter Maximum iterations
//' @param tol Convergence tolerance
//' @param verbose Print progress
//' @param n_threads Number of OpenMP threads (0 = auto)
//' @param seed Random seed (0 = random)
//' @return List with S, A, W, K, center, iterations, converged
//'
//' @keywords internal
// [[Rcpp::export]]
List cpp_fastICA_dense(const Eigen::MatrixXd& X,
                       int n_comp,
                       int alg_type,
                       int fun_type,
                       double alpha,
                       int whiten_method,
                       int max_iter,
                       double tol,
                       bool verbose,
                       int n_threads,
                       unsigned int seed) {

    // Dispatch based on nonlinearity type
    switch (fun_type) {
        case 0: {
            // LogCosh
            LogCosh<double> nonlin(alpha);
            return runICA(X, n_comp, alg_type, whiten_method, max_iter,
                          tol, verbose, n_threads, seed, nonlin);
        }
        case 1: {
            // Exp
            Exp<double> nonlin;
            return runICA(X, n_comp, alg_type, whiten_method, max_iter,
                          tol, verbose, n_threads, seed, nonlin);
        }
        case 2: {
            // Cube
            Cube<double> nonlin;
            return runICA(X, n_comp, alg_type, whiten_method, max_iter,
                          tol, verbose, n_threads, seed, nonlin);
        }
        default: {
            stop("Unknown nonlinearity type");
        }
    }

    // Should never reach here
    return List::create();
}

//' Whitening only (for testing/debugging)
//'
//' @param X Numeric matrix
//' @param n_comp Number of components
//' @param method 0 = SVD, 1 = eigen
//' @return List with X_whitened, K, K_inv, mean
//'
//' @keywords internal
// [[Rcpp::export]]
List cpp_whiten(const Eigen::MatrixXd& X, int n_comp, int method) {
    WhiteningMethod wmethod = static_cast<WhiteningMethod>(method);
    auto result = whiten<double>(X, n_comp, wmethod);

    return List::create(
        Named("X_whitened") = result.X_whitened,
        Named("K") = result.K,
        Named("K_inv") = result.K_inv,
        Named("mean") = result.mean
    );
}
