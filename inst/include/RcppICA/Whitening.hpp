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
template<typename Scalar = double>
class SVDWhitener {
public:
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    WhiteningResult<Scalar> whiten(const Matrix& X, int n_components) {
        WhiteningResult<Scalar> result;
        const int n = X.rows();  // observations
        const int m = X.cols();  // variables

        // Step 1: Center the data
        result.mean = X.colwise().mean();
        Matrix X_centered = X.rowwise() - result.mean.transpose();

        // Step 2: Compute SVD of centered data
        // X = U * S * V^T
        // Use BDCSVD for large matrices (faster than JacobiSVD)
        Eigen::BDCSVD<Matrix> svd(X_centered, Eigen::ComputeThinU | Eigen::ComputeThinV);

        // Step 3: Select top n_components
        int k = std::min(n_components, static_cast<int>(svd.singularValues().size()));
        k = std::min(k, std::min(n - 1, m));  // Can't have more components than min(n-1, m)

        // V_k: m x k (principal directions in variable space)
        Matrix V_k = svd.matrixV().leftCols(k);

        // S_k: singular values
        Vector S_k = svd.singularValues().head(k);

        // Step 4: Compute whitening matrix
        // For whitening, we want Cov(X_whitened) = I
        // X_whitened = X_centered * K where K = V * S^{-1} / sqrt(n-1)
        // This scales each PC to have unit variance
        Scalar scale_factor = std::sqrt(static_cast<Scalar>(n - 1));
        Vector scale = S_k.array().inverse() * scale_factor;
        result.K = V_k * scale.asDiagonal();

        // Dewhitening matrix: K^{-1} = V * S / sqrt(n-1)
        result.K_inv = V_k * (S_k / scale_factor).asDiagonal();

        // Step 5: Apply whitening transformation
        result.X_whitened.noalias() = X_centered * result.K;

        return result;
    }
};

// Eigendecomposition-based whitening (alternative, more memory efficient when n >> m)
template<typename Scalar = double>
class EigenWhitener {
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

        // Step 5: Apply whitening
        result.X_whitened.noalias() = X_centered * result.K;

        return result;
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
