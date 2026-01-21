#ifndef RCPPICA_ICAPARALLEL_HPP
#define RCPPICA_ICAPARALLEL_HPP

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <random>
#include "Whitening.hpp"
#include "Nonlinearities.hpp"
#include "ICADeflation.hpp"  // For ICAResult

#ifdef _OPENMP
#include <omp.h>
#endif

namespace RcppICA {

// Parallel (Symmetric) FastICA algorithm
// Extracts all independent components simultaneously with symmetric orthogonalization
template<typename Scalar = double, typename Nonlinearity = LogCosh<Scalar>>
class ICAParallel {
public:
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using RowVector = Eigen::Matrix<Scalar, 1, Eigen::Dynamic>;

private:
    int n_components_;
    int max_iter_;
    Scalar tol_;
    bool verbose_;
    int n_threads_;
    Nonlinearity nonlinearity_;
    WhiteningMethod whiten_method_;

    // Symmetric orthogonalization: W = (W * W^T)^{-1/2} * W
    void symmetricOrthogonalize(Matrix& W) {
        // Compute W * W^T (n_comp x n_comp)
        Matrix WWT;
        WWT.noalias() = W * W.transpose();

        // Eigendecomposition: WWT = V * D * V^T
        Eigen::SelfAdjointEigenSolver<Matrix> eig(WWT);

        // (W * W^T)^{-1/2} = V * D^{-1/2} * V^T
        Vector D_inv_sqrt = eig.eigenvalues().array().rsqrt();
        Matrix WWT_inv_sqrt = eig.eigenvectors()
            * D_inv_sqrt.asDiagonal()
            * eig.eigenvectors().transpose();

        // W_new = (W * W^T)^{-1/2} * W
        Matrix W_new;
        W_new.noalias() = WWT_inv_sqrt * W;
        W = W_new;
    }

    // Check convergence by comparing W matrices
    bool checkConvergence(const Matrix& W_old, const Matrix& W_new) {
        // Compare |W_new * W_old^T| diagonal to 1
        // This accounts for sign ambiguity
        Matrix product;
        product.noalias() = W_new * W_old.transpose();

        Scalar max_diff = 0;
        for (int i = 0; i < product.rows(); ++i) {
            Scalar diff = std::abs(std::abs(product(i, i)) - 1.0);
            if (diff > max_diff) max_diff = diff;
        }

        return max_diff < tol_;
    }

public:
    ICAParallel(int n_components,
                int max_iter = 200,
                Scalar tol = 1e-4,
                bool verbose = false,
                int n_threads = 0,
                Nonlinearity nonlin = Nonlinearity(),
                WhiteningMethod whiten_method = WhiteningMethod::SVD)
        : n_components_(n_components)
        , max_iter_(max_iter)
        , tol_(tol)
        , verbose_(verbose)
        , n_threads_(n_threads)
        , nonlinearity_(nonlin)
        , whiten_method_(whiten_method)
    {
#ifdef _OPENMP
        if (n_threads_ <= 0) {
            n_threads_ = omp_get_max_threads();
        }
        omp_set_num_threads(n_threads_);
#endif
    }

    ICAResult<Scalar> fit(const Matrix& X, unsigned int seed = 0) {
        const int n = X.rows();  // observations
        const int m = X.cols();  // variables

        ICAResult<Scalar> result;
        result.converged = false;
        result.iterations = 0;

        // Step 1: Whiten the data
        auto whitening = whiten<Scalar>(X, n_components_, whiten_method_);
        Matrix X_white = whitening.X_whitened;  // n x k
        result.K = whitening.K;
        result.center = whitening.mean;

        const int k = X_white.cols();  // Whitened dimension

        // Step 2: Initialize W with random orthogonal matrix
        std::mt19937 rng(seed == 0 ? std::random_device{}() : seed);
        std::normal_distribution<Scalar> dist(0.0, 1.0);

        Matrix W(n_components_, k);
        for (int i = 0; i < n_components_; ++i) {
            for (int j = 0; j < k; ++j) {
                W(i, j) = dist(rng);
            }
        }

        // Orthogonalize initial W
        symmetricOrthogonalize(W);

        // Pre-allocate matrices for the iteration
        Matrix W_new(n_components_, k);
        Matrix WX(n_components_, n);      // W * X_white^T
        Matrix GWX(n_components_, n);     // g(W * X_white^T)
        Matrix GPWX(n_components_, n);    // g'(W * X_white^T)

        // Step 3: Symmetric fixed-point iteration
        for (int iter = 0; iter < max_iter_; ++iter) {
            Matrix W_old = W;

            // Compute WX = W * X_white^T (n_comp x n)
            WX.noalias() = W * X_white.transpose();

            // Update all components (can be parallelized)
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic) num_threads(n_threads_)
#endif
            for (int p = 0; p < n_components_; ++p) {
                // Extract row p of WX
                RowVector wx_p = WX.row(p);

                // Apply nonlinearity to this component
                Matrix gwx_p, gpwx_p;
                nonlinearity_.compute(wx_p.transpose(), gwx_p, gpwx_p);

                // Store for this component
                GWX.row(p) = gwx_p.transpose();
                GPWX.row(p) = gpwx_p.transpose();
            }

            // Compute updates for all components
#ifdef _OPENMP
            #pragma omp parallel for schedule(static) num_threads(n_threads_)
#endif
            for (int p = 0; p < n_components_; ++p) {
                // Fixed-point update for component p:
                // w_new = E{X * g(w^T * X)} - E{g'(w^T * X)} * w
                // w_new = X_white^T * g(wx_p) / n - mean(g'(wx_p)) * w

                // First term: X_white^T * g(wx_p) / n
                // X_white: n x k, GWX.row(p)^T: n x 1
                // Result: k x 1
                Vector term1 = X_white.transpose() * GWX.row(p).transpose() / static_cast<Scalar>(n);

                // Second term: mean(g'(wx_p)) * w
                Scalar mean_gpwx = GPWX.row(p).mean();
                Vector w_old_p = W.row(p).transpose();
                Vector term2 = mean_gpwx * w_old_p;

                // Update
                W_new.row(p) = (term1 - term2).transpose();
            }

            // Symmetric orthogonalization
            symmetricOrthogonalize(W_new);

            W = W_new;
            ++result.iterations;

            // Check convergence
            if (checkConvergence(W_old, W)) {
                result.converged = true;
                break;
            }
        }

        if (!result.converged && verbose_) {
            Rcpp::Rcout << "Algorithm did not converge after "
                        << max_iter_ << " iterations\n";
        }

        // Store final W
        result.W = W;

        // Compute output matrices
        // S = W * X_white^T (n_components x n)
        result.S.noalias() = W * X_white.transpose();

        // A = K * W^T (mixing matrix: m x n_components)
        result.A.noalias() = result.K * W.transpose();

        return result;
    }
};

} // namespace RcppICA

#endif // RCPPICA_ICAPARALLEL_HPP
