#ifndef RCPPICA_WHITENING_HPP
#define RCPPICA_WHITENING_HPP

#include <Eigen/Dense>
#include <Eigen/SVD>

namespace RcppICA {

// Whitening method enumeration
enum class WhiteningMethod { SVD = 0, EIGEN = 1, SPECTRA = 2 };

// Result structure for whitening
template<typename Scalar = double>
struct WhiteningResult {
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> X_whitened;  // Whitened data
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> K;           // Whitening matrix
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> K_inv;       // Dewhitening matrix
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> mean;                     // Column means
};

// SVD-based whitening (recommended for numerical stability)
// Projects data onto principal components and scales to unit variance
// For large matrices, uses covariance-based approach to save memory
template<typename Scalar = double>
class SVDWhitener {
public:
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    WhiteningResult<Scalar> whiten(const Matrix& X, int n_components) {
        WhiteningResult<Scalar> result;
        const int n = X.rows();  // observations
        const int m = X.cols();  // variables

        // Step 1: Compute column means
        result.mean = X.colwise().mean();

        // For large matrices, use covariance-based approach (memory efficient)
        // Threshold: if n*m > 10 million elements, switch to covariance method
        // This avoids allocating a 40GB+ centered matrix for large datasets
        const size_t large_matrix_threshold = 10000000;
        if (static_cast<size_t>(n) * m > large_matrix_threshold) {
            return whitenLargeMatrix(X, n_components, n, m, result.mean);
        }

        // For small matrices, use direct SVD (more numerically stable)
        Matrix X_centered = X.rowwise() - result.mean.transpose();

        // Compute SVD of centered data: X = U * S * V^T
        // Use BDCSVD for large matrices (faster than JacobiSVD)
        Eigen::BDCSVD<Matrix> svd(X_centered, Eigen::ComputeThinU | Eigen::ComputeThinV);

        // Select top n_components
        int k = std::min(n_components, static_cast<int>(svd.singularValues().size()));
        k = std::min(k, std::min(n - 1, m));  // Can't have more components than min(n-1, m)

        // V_k: m x k (principal directions in variable space)
        Matrix V_k = svd.matrixV().leftCols(k);

        // S_k: singular values
        Vector S_k = svd.singularValues().head(k);

        // Compute whitening matrix
        // For whitening, we want Cov(X_whitened) = I
        // X_whitened = X_centered * K where K = V * S^{-1} / sqrt(n-1)
        Scalar scale_factor = std::sqrt(static_cast<Scalar>(n - 1));
        Vector scale = S_k.array().inverse() * scale_factor;
        result.K = V_k * scale.asDiagonal();

        // Dewhitening matrix: K^{-1} = V * S / sqrt(n-1)
        result.K_inv = V_k * (S_k / scale_factor).asDiagonal();

        // Apply whitening transformation
        result.X_whitened.noalias() = X_centered * result.K;

        return result;
    }

private:
    // Memory-efficient whitening for large matrices using covariance approach
    WhiteningResult<Scalar> whitenLargeMatrix(
        const Matrix& X,
        int n_components,
        int n, int m,
        const Vector& mean)
    {
        WhiteningResult<Scalar> result;
        result.mean = mean;

        // Compute covariance WITHOUT materializing X_centered
        // Identity: Cov(X) = [X'X - n*μμ'] / (n-1)
        Matrix XtX;
        XtX.noalias() = X.transpose() * X;  // m×m, delegates to BLAS

        Matrix mean_outer;
        mean_outer.noalias() = mean * mean.transpose();

        const Scalar n_scalar = static_cast<Scalar>(n);
        const Scalar nm1 = static_cast<Scalar>(n - 1);
        Matrix C = XtX / nm1;
        C -= (n_scalar / nm1) * mean_outer;

        // Eigendecomposition of covariance (equivalent to SVD for whitening)
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

        // Whitening matrix: K = E * D^{-1/2}
        Vector scale = D_k.array().rsqrt();
        result.K = E_k * scale.asDiagonal();

        // Dewhitening: K^{-1} = E * D^{1/2}
        result.K_inv = E_k * D_k.array().sqrt().matrix().asDiagonal();

        // Apply whitening in chunks
        result.X_whitened = applyWhiteningChunked(X, mean, result.K, n, m, k);

        return result;
    }

    // Apply whitening in chunks to avoid full matrix allocation
    Matrix applyWhiteningChunked(
        const Matrix& X,
        const Vector& mean,
        const Matrix& K,
        int n, int m, int k)
    {
        Matrix X_whitened(n, k);

        // Target: ~100MB per chunk
        const int target_chunk_mb = 100;
        const int chunk_size = std::max(1000, std::min(
            static_cast<int>((target_chunk_mb * 1024 * 1024) / (8 * m)),
            n
        ));

        for (int start = 0; start < n; start += chunk_size) {
            int end = std::min(start + chunk_size, n);
            int chunk_rows = end - start;

            Matrix X_chunk_centered = X.middleRows(start, chunk_rows).rowwise() - mean.transpose();
            X_whitened.middleRows(start, chunk_rows).noalias() = X_chunk_centered * K;
        }

        return X_whitened;
    }
};

// Eigendecomposition-based whitening (alternative, more memory efficient when n >> m)
// Optimized to avoid materializing the full n×m centered matrix
template<typename Scalar = double>
class EigenWhitener {
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

        // Step 3: Eigendecomposition of symmetric covariance matrix
        Eigen::SelfAdjointEigenSolver<Matrix> eig(C);

        // Eigenvalues are in ascending order, we want descending (largest first)
        int k = std::min(n_components, m);

        // Take top k eigenvectors (from end, since ascending order)
        // E_k: m x k (eigenvectors as columns)
        Matrix E_k = eig.eigenvectors().rightCols(k).rowwise().reverse();
        Vector D_k = eig.eigenvalues().tail(k).reverse();

        // Handle small/negative eigenvalues (numerical issues)
        const Scalar eps = std::numeric_limits<Scalar>::epsilon() * 100;
        for (int i = 0; i < k; ++i) {
            if (D_k(i) < eps) {
                D_k(i) = eps;
            }
        }

        // Step 4: Whitening matrix: K = E * D^{-1/2}
        Vector scale = D_k.array().rsqrt();
        result.K = E_k * scale.asDiagonal();

        // Dewhitening: K^{-1} = E * D^{1/2}
        result.K_inv = E_k * D_k.array().sqrt().matrix().asDiagonal();

        // Step 5: Apply whitening in chunks to avoid full X_centered allocation
        result.X_whitened = applyWhiteningChunked(X, result.mean, result.K, n, m, k);

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
};

} // namespace RcppICA

// Include Spectra and sparse whiteners after namespace declaration
#include "WhiteningSpectra.hpp"
#include "WhiteningSparseCov.hpp"

namespace RcppICA {

// Factory function for whitening (dense matrices)
template<typename Scalar = double>
WhiteningResult<Scalar> whiten(
    const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& X,
    int n_components,
    WhiteningMethod method = WhiteningMethod::SPECTRA)  // Changed default to SPECTRA
{
    if (method == WhiteningMethod::SPECTRA) {
        // SPECTRA method (default, fastest - computes only top-k eigenvalues)
        SpectraWhitener<Scalar> whitener;
        return whitener.whiten(X, n_components);
    } else if (method == WhiteningMethod::SVD) {
        SVDWhitener<Scalar> whitener;
        return whitener.whiten(X, n_components);
    } else {
        // EIGEN method (full eigendecomposition)
        EigenWhitener<Scalar> whitener;
        return whitener.whiten(X, n_components);
    }
}

// Factory function for whitening (sparse matrices)
// For sparse input, always use sparse-aware covariance computation
// This avoids materializing the dense centered matrix, saving massive memory
template<typename Scalar = double>
WhiteningResult<Scalar> whiten(
    const Eigen::SparseMatrix<Scalar>& X_sparse,
    int n_components,
    WhiteningMethod method = WhiteningMethod::SPECTRA)
{
    // For sparse matrices, use sparse-aware covariance computation
    // Ignore method parameter - sparse covariance is always optimal
    SparseWhitener<Scalar> whitener;
    return whitener.whiten(X_sparse, n_components);
}

} // namespace RcppICA

#endif // RCPPICA_WHITENING_HPP
