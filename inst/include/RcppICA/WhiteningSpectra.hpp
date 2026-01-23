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
// Optimized to avoid materializing the full n×m centered matrix
template<typename Scalar = double>
class SpectraWhitener {
public:
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    WhiteningResult<Scalar> whiten(const Matrix& X, int n_components) {
        WhiteningResult<Scalar> result;
        const int n = X.rows();
        const int m = X.cols();

        // Step 1: Compute column means
        result.mean = X.colwise().mean();

        // Step 2: Compute covariance WITHOUT materializing X_centered
        // Mathematical identity: Cov(X) = [X'X - n*μμ'] / (n-1)
        // This avoids allocating the n×m centered matrix (massive memory savings)
        Matrix XtX;
        XtX.noalias() = X.transpose() * X;  // m×m, delegates to BLAS DGEMM

        // Compute mean outer product: μμ' (m×m, rank-1)
        Matrix mean_outer;
        mean_outer.noalias() = result.mean * result.mean.transpose();

        // Combine: C = XtX/(n-1) - n/(n-1) * μμ'
        const Scalar n_scalar = static_cast<Scalar>(n);
        const Scalar nm1 = static_cast<Scalar>(n - 1);
        Matrix C = XtX / nm1;
        C -= (n_scalar / nm1) * mean_outer;

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
                return fallbackToEigen(X, result.mean, C, n_components, n, m);
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

            // Step 6: Apply whitening in chunks to avoid full X_centered allocation
            result.X_whitened = applyWhiteningChunked(X, result.mean, result.K, n, m, k);

        } catch (const std::exception& e) {
            // Fallback to Eigen if Spectra throws any exception
            return fallbackToEigen(X, result.mean, C, n_components, n, m);
        }

        return result;
    }

private:
    // Apply whitening in chunks to avoid full matrix densification
    Matrix applyWhiteningChunked(
        const Matrix& X,
        const Vector& mean,
        const Matrix& K,
        int n, int m, int k)
    {
        Matrix X_whitened(n, k);

        // Determine chunk size based on available memory
        // Target: ~100MB per chunk at double precision
        // 100MB / 8 bytes / m columns ≈ chunk_size rows
        const int target_chunk_mb = 100;
        const int chunk_size = std::max(1000, std::min(
            static_cast<int>((target_chunk_mb * 1024 * 1024) / (8 * m)),
            n
        ));

        // Process data in chunks
        for (int start = 0; start < n; start += chunk_size) {
            int end = std::min(start + chunk_size, n);
            int chunk_rows = end - start;

            // Center chunk and apply whitening
            Matrix X_chunk_centered = X.middleRows(start, chunk_rows).rowwise() - mean.transpose();

            // Apply whitening to chunk
            X_whitened.middleRows(start, chunk_rows).noalias() = X_chunk_centered * K;
        }

        return X_whitened;
    }

    // Fallback to full eigendecomposition if Spectra fails
    WhiteningResult<Scalar> fallbackToEigen(
        const Matrix& X,
        const Vector& mean,
        const Matrix& C,
        int n_components,
        int n,
        int m)
    {
        WhiteningResult<Scalar> result;
        result.mean = mean;

        // Full eigendecomposition (covariance C already computed)
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

        // Apply whitening in chunks
        result.X_whitened = applyWhiteningChunked(X, mean, result.K, n, m, k);

        return result;
    }
};

} // namespace RcppICA

#endif // RCPPICA_WHITENING_SPECTRA_HPP
