#ifndef RCPPICA_ICADEFLATION_HPP
#define RCPPICA_ICADEFLATION_HPP

#include <Eigen/Dense>
#include <random>
#include "Whitening.hpp"
#include "Nonlinearities.hpp"

namespace RcppICA {

// ICA Result structure
template<typename Scalar = double>
struct ICAResult {
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> S;       // Independent components (n_comp x n)
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> A;       // Mixing matrix (m x n_comp)
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> W;       // Unmixing matrix (n_comp x k)
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> K;       // Whitening matrix (m x k)
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> center;               // Column means
    int iterations;
    bool converged;
};

// Deflation-based FastICA algorithm
// Extracts independent components one at a time
template<typename Scalar = double, typename Nonlinearity = LogCosh<Scalar>>
class ICADeflation {
public:
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using RowVector = Eigen::Matrix<Scalar, 1, Eigen::Dynamic>;

private:
    int n_components_;
    int max_iter_;
    Scalar tol_;
    bool verbose_;
    Nonlinearity nonlinearity_;
    WhiteningMethod whiten_method_;

public:
    ICADeflation(int n_components,
                 int max_iter = 200,
                 Scalar tol = 1e-4,
                 bool verbose = false,
                 Nonlinearity nonlin = Nonlinearity(),
                 WhiteningMethod whiten_method = WhiteningMethod::SVD)
        : n_components_(n_components)
        , max_iter_(max_iter)
        , tol_(tol)
        , verbose_(verbose)
        , nonlinearity_(nonlin)
        , whiten_method_(whiten_method)
    {}

    ICAResult<Scalar> fit(const Matrix& X, unsigned int seed = 0) {
        const int n = X.rows();  // observations
        const int m = X.cols();  // variables

        ICAResult<Scalar> result;
        result.converged = true;
        result.iterations = 0;

        // Step 1: Whiten the data
        auto whitening = whiten<Scalar>(X, n_components_, whiten_method_);
        Matrix X_white = whitening.X_whitened;  // n x k
        result.K = whitening.K;
        result.center = whitening.mean;

        const int k = X_white.cols();  // Whitened dimension

        // Step 2: Initialize random number generator
        std::mt19937 rng(seed == 0 ? std::random_device{}() : seed);
        std::normal_distribution<Scalar> dist(0.0, 1.0);

        // Step 3: Initialize W matrix (n_components x k)
        result.W.resize(n_components_, k);

        // Step 4: Extract components one by one (deflation)
        for (int p = 0; p < n_components_; ++p) {
            // Random initialization for this component
            Vector w(k);
            for (int i = 0; i < k; ++i) {
                w(i) = dist(rng);
            }
            w.normalize();

            int iter = 0;
            bool converged = false;

            while (iter < max_iter_ && !converged) {
                Vector w_old = w;

                // Fixed-point iteration:
                // w_new = E{X * g(w^T * X)} - E{g'(w^T * X)} * w

                // Compute w^T * X for all observations
                // X_white: n x k, w: k x 1
                // wx: n x 1
                Vector wx = X_white * w;

                // Apply nonlinearity g(wx) and g'(wx)
                Matrix gwx, gpwx;
                nonlinearity_.compute(wx, gwx, gpwx);

                // First term: E{X * g(w^T * X)} = X^T * g(wx) / n
                // X_white^T: k x n, gwx: n x 1
                // Result: k x 1
                Vector term1 = X_white.transpose() * gwx / static_cast<Scalar>(n);

                // Second term: E{g'(w^T * X)} * w = mean(g'(wx)) * w
                Scalar mean_gpwx = gpwx.mean();
                Vector term2 = mean_gpwx * w;

                // Update
                w = term1 - term2;

                // Gram-Schmidt orthogonalization against previous components
                for (int j = 0; j < p; ++j) {
                    Vector wj = result.W.row(j).transpose();
                    w -= w.dot(wj) * wj;
                }

                // Normalize
                Scalar norm = w.norm();
                if (norm < std::numeric_limits<Scalar>::epsilon() * 100) {
                    // Reinitialize if norm is too small
                    for (int i = 0; i < k; ++i) {
                        w(i) = dist(rng);
                    }
                    w.normalize();
                } else {
                    w /= norm;
                }

                // Check convergence (allowing for sign flip)
                Scalar change = std::min((w - w_old).norm(), (w + w_old).norm());
                converged = (change < tol_);

                ++iter;
            }

            result.W.row(p) = w.transpose();
            result.iterations = std::max(result.iterations, iter);

            if (!converged) {
                result.converged = false;
                if (verbose_) {
                    Rcpp::Rcout << "Component " << p + 1
                                << " did not converge after "
                                << max_iter_ << " iterations\n";
                }
            }
        }

        // Compute output matrices
        // S = W * X_white^T (n_components x n)
        result.S.noalias() = result.W * X_white.transpose();

        // A = K * W^T (mixing matrix: m x n_components)
        // This transforms from IC space back to original variable space
        result.A.noalias() = result.K * result.W.transpose();

        return result;
    }
};

} // namespace RcppICA

#endif // RCPPICA_ICADEFLATION_HPP
