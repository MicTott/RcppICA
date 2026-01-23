#ifndef RCPPICA_WHITENING_SPARSE_COV_HPP
#define RCPPICA_WHITENING_SPARSE_COV_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/Util/SelectionRule.h>
#include <Spectra/Util/CompInfo.h>

namespace RcppICA {

// Forward declaration of WhiteningResult
template<typename Scalar>
struct WhiteningResult;

// Sparse-aware whitening using streaming covariance computation
// Computes X'X without materializing centered matrix (major memory savings)
template<typename Scalar = double>
class SparseWhitener {
public:
    using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    using SparseMatrix = Eigen::SparseMatrix<Scalar>;

    WhiteningResult<Scalar> whiten(const SparseMatrix& X_sparse, int n_components) {
        WhiteningResult<Scalar> result;
        const int n = X_sparse.rows();
        const int m = X_sparse.cols();

        // Step 1: Compute column means (sparse-aware iteration)
        result.mean = Vector::Zero(m);

        // Iterate over non-zero entries only
        for (int k = 0; k < X_sparse.outerSize(); ++k) {
            for (typename SparseMatrix::InnerIterator it(X_sparse, k); it; ++it) {
                result.mean(it.col()) += it.value();
            }
        }
        result.mean /= static_cast<Scalar>(n);

        // Step 2: Compute X'X
        // For highly sparse data, sparse multiplication is faster.
        // For moderate sparsity (>5% non-zero), BLAS dense multiply wins
        // due to vectorization and cache-optimal memory access.
        double density = static_cast<double>(X_sparse.nonZeros()) /
                         (static_cast<double>(n) * m);

        Matrix XtX_raw(m, m);
        if (density < 0.05) {
            // Sparse path: only beneficial at high sparsity
            SparseMatrix Xt = X_sparse.transpose();
            XtX_raw = Matrix(Xt * X_sparse);
        } else {
            // Dense path: convert and use BLAS DGEMM (faster for moderate sparsity)
            Matrix X_dense = Matrix(X_sparse);
            XtX_raw.noalias() = X_dense.transpose() * X_dense;
        }

        // Step 3: Compute centered covariance without materializing X_centered
        // Mathematical insight: C = (X - μ)'(X - μ) / (n-1)
        //                          = [X'X - n*μμ'] / (n-1)
        // This avoids creating the dense n×m centered matrix
        Matrix mean_outer = result.mean * result.mean.transpose();
        const Scalar n_scalar = static_cast<Scalar>(n);
        const Scalar nm1 = static_cast<Scalar>(n - 1);
        Matrix C = XtX_raw / nm1;
        C -= (n_scalar / nm1) * mean_outer;

        // Step 4: Use Spectra for top-k eigenvalues (unchanged from dense version)
        int k = std::min(n_components, m);
        int ncv = std::min(m, std::max(2 * k + 1, std::min(20, m)));

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
                // Fallback to full eigendecomposition if Spectra fails
                return fallbackToEigen(X_sparse, result.mean, C, n_components, n, m);
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

            // Step 6: Apply whitening using chunked processing to avoid full densification
            // For very large n, process in chunks to limit memory usage
            result.X_whitened = applyWhiteningChunked(X_sparse, result.mean, result.K, n, m, k);

        } catch (const std::exception& e) {
            // Fallback to Eigen if Spectra throws any exception
            return fallbackToEigen(X_sparse, result.mean, C, n_components, n, m);
        }

        return result;
    }

private:
    // Apply whitening in chunks to avoid full matrix densification
    Matrix applyWhiteningChunked(
        const SparseMatrix& X_sparse,
        const Vector& mean,
        const Matrix& K,
        int n, int m, int k)
    {
        Matrix X_whitened(n, k);

        // Precompute mean contribution: mean' * K (1 x k vector)
        // This lets us avoid centering each chunk explicitly:
        //   (X_chunk - mean) * K = X_chunk * K - mean' * K
        Eigen::Matrix<Scalar, 1, Eigen::Dynamic> meanK = mean.transpose() * K;

        // Determine chunk size based on available memory
        // Target: ~100MB per chunk at double precision
        const int target_chunk_mb = 100;
        const int chunk_size = std::max(1000, std::min(
            static_cast<int>((target_chunk_mb * 1024 * 1024) / (8 * m)),
            n
        ));

        // Process data in chunks
        for (int start = 0; start < n; start += chunk_size) {
            int end = std::min(start + chunk_size, n);
            int chunk_rows = end - start;

            // Convert sparse chunk to dense and apply whitening in one step
            // (X_chunk - mean) * K = X_chunk * K - meanK
            Matrix X_chunk = Matrix(X_sparse.middleRows(start, chunk_rows));
            X_whitened.middleRows(start, chunk_rows).noalias() = X_chunk * K;
            X_whitened.middleRows(start, chunk_rows).rowwise() -= meanK;
        }

        return X_whitened;
    }

    // Fallback to full eigendecomposition if Spectra fails
    WhiteningResult<Scalar> fallbackToEigen(
        const SparseMatrix& X_sparse,
        const Vector& mean,
        const Matrix& C,
        int n_components,
        int n,
        int m)
    {
        WhiteningResult<Scalar> result;
        result.mean = mean;

        // Full eigendecomposition of covariance
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

        // Apply whitening (chunked)
        result.X_whitened = applyWhiteningChunked(X_sparse, mean, result.K, n, m, k);

        return result;
    }
};

} // namespace RcppICA

#endif // RCPPICA_WHITENING_SPARSE_COV_HPP
