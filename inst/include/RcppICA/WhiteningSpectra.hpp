#ifndef RCPPICA_WHITENING_SPECTRA_HPP
#define RCPPICA_WHITENING_SPECTRA_HPP

#include <Eigen/Dense>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/Util/SelectionRule.h>
#include <Spectra/Util/CompInfo.h>

namespace RcppICA {

// Forward declaration of WhiteningResult
template<typename Scalar>
struct WhiteningResult;

// Spectra-based whitening (fastest - computes only top-k eigenvalues)
// Uses Lanczos iteration for partial eigendecomposition
template<typename Scalar = double>
class SpectraWhitener {
public:
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    WhiteningResult<Scalar> whiten(const Matrix& X, int n_components) {
        WhiteningResult<Scalar> result;
        const int n = X.rows();
        const int m = X.cols();

        // Step 1: Center the data
        result.mean = X.colwise().mean();
        Matrix X_centered = X.rowwise() - result.mean.transpose();

        // Step 2: Compute covariance matrix: C = X^T * X / (n-1)
        // This is m x m
        Matrix C;
        C.noalias() = X_centered.transpose() * X_centered;
        C /= static_cast<Scalar>(n - 1);

        // Step 3: Determine number of eigenvalues to compute
        int k = std::min(n_components, m);

        // For Spectra, we need ncv (number of Lanczos basis vectors)
        // ncv must satisfy: k < ncv <= m
        // Recommended: ncv >= 2*k for good convergence
        int ncv = std::min(m, std::max(2 * k + 1, std::min(20, m)));

        // Step 4: Use Spectra for partial eigendecomposition
        // Only compute top-k eigenvalues instead of all m
        try {
            // Create matrix operation object for Spectra
            Spectra::DenseSymMatProd<Scalar> op(C);

            // Create eigen solver for k largest eigenvalues
            Spectra::SymEigsSolver<Spectra::DenseSymMatProd<Scalar>> eigs(op, k, ncv);

            // Initialize with random vector and compute
            eigs.init();
            eigs.compute(Spectra::SortRule::LargestAlge, 1000, 1e-10);

            // Check convergence
            if (eigs.info() != Spectra::CompInfo::Successful) {
                // Fallback to Eigen if Spectra fails to converge
                return fallbackToEigen(X, X_centered, n_components, n, m);
            }

            // Extract eigenvalues and eigenvectors (already in descending order)
            Vector eigenvalues = eigs.eigenvalues();
            Matrix eigenvectors = eigs.eigenvectors();

            // Handle small/negative eigenvalues (numerical issues)
            const Scalar eps = std::numeric_limits<Scalar>::epsilon() * 100;
            for (int i = 0; i < k; ++i) {
                if (eigenvalues(i) < eps) {
                    eigenvalues(i) = eps;
                }
            }

            // Step 5: Whitening matrix: K = E * D^{-1/2}
            Vector scale = eigenvalues.array().rsqrt();
            result.K = eigenvectors * scale.asDiagonal();

            // Dewhitening: K^{-1} = E * D^{1/2}
            result.K_inv = eigenvectors * eigenvalues.array().sqrt().matrix().asDiagonal();

            // Step 6: Apply whitening
            result.X_whitened.noalias() = X_centered * result.K;

        } catch (const std::exception& e) {
            // Fallback to Eigen if Spectra throws any exception
            return fallbackToEigen(X, X_centered, n_components, n, m);
        }

        return result;
    }

private:
    // Fallback to full eigendecomposition if Spectra fails
    WhiteningResult<Scalar> fallbackToEigen(
        const Matrix& X,
        const Matrix& X_centered,
        int n_components,
        int n,
        int m)
    {
        WhiteningResult<Scalar> result;
        result.mean = X.colwise().mean();

        // Compute covariance matrix
        Matrix C;
        C.noalias() = X_centered.transpose() * X_centered;
        C /= static_cast<Scalar>(n - 1);

        // Full eigendecomposition
        Eigen::SelfAdjointEigenSolver<Matrix> eig(C);

        int k = std::min(n_components, m);

        // Take top k eigenvectors (from end, since ascending order)
        Matrix E_k = eig.eigenvectors().rightCols(k).rowwise().reverse();
        Vector D_k = eig.eigenvalues().tail(k).reverse();

        // Handle small/negative eigenvalues
        const Scalar eps = std::numeric_limits<Scalar>::epsilon() * 100;
        for (int i = 0; i < k; ++i) {
            if (D_k(i) < eps) {
                D_k(i) = eps;
            }
        }

        // Whitening matrix
        Vector scale = D_k.array().rsqrt();
        result.K = E_k * scale.asDiagonal();

        // Dewhitening
        result.K_inv = E_k * D_k.array().sqrt().matrix().asDiagonal();

        // Apply whitening
        result.X_whitened.noalias() = X_centered * result.K;

        return result;
    }
};

} // namespace RcppICA

#endif // RCPPICA_WHITENING_SPECTRA_HPP
