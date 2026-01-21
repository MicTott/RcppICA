// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

#include "../inst/include/RcppICA.h"
#include <RcppEigen.h>

using namespace Rcpp;
using namespace RcppICA;

// Helper function to run ICA with a specific nonlinearity (sparse version)
// The key difference: takes sparse matrix, whitens using sparse-aware method,
// then runs standard ICA on the (dense) whitened data
template<typename Nonlinearity>
List runICASparse(const Eigen::SparseMatrix<double>& X_sparse,
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

    // Step 1: Whiten using sparse-aware method (major memory savings here!)
    auto whitening = whiten<double>(X_sparse, n_comp, wmethod);

    // From here on, we work with dense whitened data (which is small: n x k)
    Eigen::MatrixXd X_white = whitening.X_whitened;

    // Step 2: Run ICA on whitened data (same as dense version)
    if (alg_type == 0) {
        // Parallel (symmetric) algorithm
        // We need to create a "whitening-only" version since data is already whitened
        // For now, we'll manually run the ICA iteration

        const int n = X_white.rows();
        const int k = X_white.cols();

        // Initialize W with random orthogonal matrix
        std::mt19937 rng(seed == 0 ? std::random_device{}() : seed);
        std::normal_distribution<double> dist(0.0, 1.0);

        Eigen::MatrixXd W(n_comp, k);
        for (int i = 0; i < n_comp; ++i) {
            for (int j = 0; j < k; ++j) {
                W(i, j) = dist(rng);
            }
        }

        // Symmetric orthogonalization
        auto symmetricOrthogonalize = [](Eigen::MatrixXd& W) {
            Eigen::MatrixXd WWT = W * W.transpose();
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(WWT);
            Eigen::VectorXd D_inv_sqrt = eig.eigenvalues().array().rsqrt();
            Eigen::MatrixXd WWT_inv_sqrt = eig.eigenvectors()
                * D_inv_sqrt.asDiagonal()
                * eig.eigenvectors().transpose();
            W = WWT_inv_sqrt * W;
        };

        symmetricOrthogonalize(W);

        // Fixed-point iteration
        bool converged = false;
        int iter;

        for (iter = 0; iter < max_iter; ++iter) {
            Eigen::MatrixXd W_old = W;

            // Compute WX = W * X_white^T
            Eigen::MatrixXd WX = W * X_white.transpose();

            // Update all components
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic) num_threads(n_threads)
#endif
            for (int p = 0; p < n_comp; ++p) {
                Eigen::VectorXd wx_p = WX.row(p);

                // Apply nonlinearity
                Eigen::MatrixXd gwx_p, gpwx_p;
                nonlin.compute(wx_p, gwx_p, gpwx_p);

                // Update rule: w_new = E[x*g(w'x)] - E[g'(w'x)]*w
                Eigen::VectorXd term1 = X_white.transpose() * gwx_p / n;
                double term2 = gpwx_p.mean();
                W.row(p) = term1 - term2 * W_old.row(p).transpose();
            }

            // Symmetric orthogonalization
            symmetricOrthogonalize(W);

            // Check convergence
            Eigen::MatrixXd product = W * W_old.transpose();
            double max_diff = 0;
            for (int i = 0; i < product.rows(); ++i) {
                double diff = std::abs(std::abs(product(i, i)) - 1.0);
                if (diff > max_diff) max_diff = diff;
            }

            if (max_diff < tol) {
                converged = true;
                break;
            }

            if (verbose && (iter % 10 == 0)) {
                Rcpp::Rcout << "Iteration " << iter << ", convergence: " << max_diff << std::endl;
            }
        }

        result.W = W;
        result.iterations = iter + 1;
        result.converged = converged;

    } else {
        // Deflation algorithm
        const int n = X_white.rows();
        const int k = X_white.cols();

        Eigen::MatrixXd W(n_comp, k);
        bool converged = true;
        int total_iter = 0;

        for (int comp = 0; comp < n_comp; ++comp) {
            // Initialize w
            std::mt19937 rng(seed == 0 ? std::random_device{}() : seed + comp);
            std::normal_distribution<double> dist(0.0, 1.0);

            Eigen::VectorXd w(k);
            for (int j = 0; j < k; ++j) {
                w(j) = dist(rng);
            }
            w.normalize();

            bool comp_converged = false;
            int iter;

            for (iter = 0; iter < max_iter; ++iter) {
                Eigen::VectorXd w_old = w;

                // Compute w'X
                Eigen::VectorXd wx = X_white * w;

                // Apply nonlinearity
                Eigen::MatrixXd gwx, gpwx;
                nonlin.compute(wx, gwx, gpwx);

                // Update: w_new = E[x*g(w'x)] - E[g'(w'x)]*w
                Eigen::VectorXd term1 = X_white.transpose() * gwx / n;
                double term2 = gpwx.mean();
                w = term1 - term2 * w_old;

                // Orthogonalize against previous components
                for (int j = 0; j < comp; ++j) {
                    w -= w.dot(W.row(j)) * W.row(j).transpose();
                }

                w.normalize();

                // Check convergence
                double diff = std::abs(std::abs(w.dot(w_old)) - 1.0);
                if (diff < tol) {
                    comp_converged = true;
                    break;
                }
            }

            W.row(comp) = w;
            total_iter += iter + 1;

            if (!comp_converged) {
                converged = false;
            }

            if (verbose) {
                Rcpp::Rcout << "Component " << comp + 1 << " converged in "
                           << iter + 1 << " iterations" << std::endl;
            }
        }

        result.W = W;
        result.iterations = total_iter;
        result.converged = converged;
    }

    // Step 3: Compute sources and mixing matrix
    result.S = result.W * X_white.transpose();  // n_comp x n
    result.A = whitening.K * result.W.transpose();  // m x n_comp
    result.K = whitening.K;
    result.center = whitening.mean;

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

//' Fast Independent Component Analysis (Sparse Matrix Implementation)
//'
//' @param X_sparse Sparse matrix (dgCMatrix, observations x variables)
//' @param n_comp Number of independent components to extract
//' @param alg_type Algorithm type: 0 = parallel (symmetric), 1 = deflation
//' @param fun_type Nonlinearity: 0 = logcosh, 1 = exp, 2 = cube
//' @param alpha Parameter for logcosh (between 1 and 2)
//' @param whiten_method Whitening method (ignored for sparse, always uses sparse covariance)
//' @param max_iter Maximum iterations
//' @param tol Convergence tolerance
//' @param verbose Print progress
//' @param n_threads Number of OpenMP threads (0 = auto)
//' @param seed Random seed (0 = random)
//' @return List with S, A, W, K, center, iterations, converged
//'
//' @keywords internal
// [[Rcpp::export]]
List cpp_fastICA_sparse(const Eigen::SparseMatrix<double>& X_sparse,
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
            return runICASparse(X_sparse, n_comp, alg_type, whiten_method, max_iter,
                               tol, verbose, n_threads, seed, nonlin);
        }
        case 1: {
            // Exp
            Exp<double> nonlin;
            return runICASparse(X_sparse, n_comp, alg_type, whiten_method, max_iter,
                               tol, verbose, n_threads, seed, nonlin);
        }
        case 2: {
            // Cube
            Cube<double> nonlin;
            return runICASparse(X_sparse, n_comp, alg_type, whiten_method, max_iter,
                               tol, verbose, n_threads, seed, nonlin);
        }
        default: {
            stop("Unknown nonlinearity type");
        }
    }

    // Should never reach here
    return List::create();
}

//' Whitening for sparse matrices (for testing/debugging)
//'
//' @param X_sparse Sparse matrix (dgCMatrix)
//' @param n_comp Number of components
//' @param method Whitening method (ignored, always uses sparse covariance)
//' @return List with X_whitened, K, K_inv, mean
//'
//' @keywords internal
// [[Rcpp::export]]
List cpp_whiten_sparse(const Eigen::SparseMatrix<double>& X_sparse,
                       int n_comp,
                       int method) {
    WhiteningMethod wmethod = static_cast<WhiteningMethod>(method);
    auto result = whiten<double>(X_sparse, n_comp, wmethod);

    return List::create(
        Named("X_whitened") = result.X_whitened,
        Named("K") = result.K,
        Named("K_inv") = result.K_inv,
        Named("mean") = result.mean
    );
}
