# RcppICA 0.2.1

## Performance

* **Eliminated dense X_centered allocation in whitening**: All three whitening
  methods (Spectra, Eigen, SVD) now compute covariance via the identity
  `Cov(X) = [X'X - n*mu*mu'] / (n-1)`, avoiding the n x m centered matrix
  allocation that dominated memory at large scales.

* **Chunked whitening application**: Whitening is applied in ~100MB chunks
  rather than allocating a full n x k result matrix upfront, bounding peak
  memory regardless of dataset size.

* **BLAS delegation enabled**: Added `EIGEN_USE_BLAS` compilation flag so
  Eigen delegates large matrix multiplications to vendor-optimized BLAS
  (Accelerate on macOS, MKL on Linux/Windows).

* **Removed EIGEN_INITIALIZE_MATRICES_BY_ZERO**: Eliminated unnecessary
  zero-initialization of all matrix allocations.

* **Sparse whitening improvements**:
  - Density-adaptive covariance: uses BLAS dense multiply when density > 5%,
    sparse multiply only at high sparsity (< 5% non-zero)
  - Eliminated double-copy in chunked whitening via precomputed mean
    contribution vector

## Documentation

* Fixed vignettes: corrected default whitening method references, removed
  inaccurate performance claims, fixed broken code examples.

* Removed vignettes with unresolvable dependencies (ica-vs-nmf-modules,
  snRNAseq-workflow).

* Cleaned up README: accurate performance table, corrected feature
  descriptions.

## Internal

* Reduced Suggests dependencies (removed igraph, scRNAseq, scran, ggplot2,
  pheatmap, patchwork that were only needed by removed vignettes).

# RcppICA 0.2.0

## New Features

* **Sparse matrix support**: Pass `dgCMatrix` objects directly to `fastICA()`.
  Covariance is computed without materializing the dense centered matrix,
  enabling ICA on large sparse datasets.

* **Spectra library integration**: Truncated eigendecomposition via Lanczos
  iteration computes only the top-k eigenvalues needed for whitening.

* **Three whitening methods**:
  - `"spectra"` (default): Truncated eigendecomposition, fastest for k << m
  - `"eigen"`: Full eigendecomposition of covariance matrix
  - `"svd"`: SVD-based, most numerically stable

# RcppICA 0.1.0

* Initial release.

* **FastICA algorithm** with parallel (symmetric) and deflation modes.

* **Nonlinearity functions**: logcosh, exp, cube.

* **OpenMP parallelization** for the parallel algorithm.

* **S4 class interface** (`ICAResult`) with print, summary, and predict methods.

* **Reproducibility** via seed parameter.
