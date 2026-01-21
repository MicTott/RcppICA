# RcppICA Performance Analysis

## Executive Summary

**RcppICA is currently 2-3x SLOWER than base fastICA**, not faster. This contradicts the original design goals and must be fixed before release.

### Benchmark Results (Fair Comparison)

| Data Size | RcppICA (eigen) | Base fastICA (C) | **Speedup** |
|-----------|-----------------|------------------|-------------|
| 500 × 100, k=5 | 0.012s | 0.009s | **0.75x (SLOWER)** |
| 1000 × 500, k=5 | 0.105s | 0.047s | **0.45x (SLOWER)** |
| 5000 × 1000, k=10 | 1.267s | 0.426s | **0.34x (SLOWER)** |

The performance gap **increases with data size**, indicating algorithmic inefficiency, not just overhead.

## Root Cause: Inefficient Whitening

The bottleneck is in the whitening step, specifically the eigendecomposition computation in [Whitening.hpp:90-94](inst/include/RcppICA/Whitening.hpp#L90-L94):

```cpp
// Compute covariance matrix: C = X^T * X / (n-1)
Matrix C;
C.noalias() = X_centered.transpose() * X_centered;
C /= static_cast<Scalar>(n - 1);

// Eigendecomposition of symmetric covariance matrix
Eigen::SelfAdjointEigenSolver<Matrix> eig(C);
```

### Why is it slow?

1. **We compute the covariance matrix explicitly**: `C = X^T * X / (n-1)` creates an m×m dense matrix
2. **For 1000 variables**: This is a 1000×1000 matrix requiring significant memory and computation
3. **SelfAdjointEigenSolver computes ALL eigenvalues/eigenvectors** (all 1000 of them)
4. **We only use the top k** (typically 5-10), discarding 99% of the computation
5. **Base fastICA uses LAPACK more efficiently**, likely with truncated methods or better-tuned algorithms

### Detailed Breakdown (1000×500 matrix, 5 components)

| Method | Time (10 calls) | Time per call | Relative |
|--------|-----------------|---------------|----------|
| **RcppICA SVD** | 1.381s | 0.138s | 3.77x slower |
| **RcppICA Eigen** | 1.026s | 0.103s | 2.81x slower |
| Base R `svd()` | 0.701s | 0.070s | 1.91x slower |
| **Base R `eigen()`** | 0.366s | 0.037s | **Baseline** |

Even base R's `eigen()` function is **2.8x faster** than our Eigen implementation, which suggests we're not using Eigen optimally.

### Why base fastICA is faster

Looking at the fastICA C code (from the CRAN package), they use:
1. **R's built-in LAPACK calls** for eigendecomposition (highly optimized, vendor-tuned BLAS)
2. **Sparse/truncated methods** that avoid computing all eigenvalues
3. **Better memory layout** for cache efficiency
4. **Optimized matrix multiplication** using BLAS level 3 operations

## Impact on Total Runtime

Whitening consumes **~95-98% of total RcppICA runtime**:

```
Total RcppICA time: 0.105s (medium dataset)
  Whitening:        0.103s (98%)
  ICA iterations:   0.002s (2%)
```

The **ICA algorithm itself is very fast** - the problem is purely in the whitening preprocessing.

## Why Our Original Claim Was Wrong

The original plan claimed "30-150x speedup" based on:
1. **Incorrect assumption**: That Eigen would be faster than R's LAPACK
2. **Focus on ICA iterations**: The iterations ARE fast, but they're not the bottleneck
3. **Missing benchmarks**: No actual comparison was performed before claiming speedup

Reality:
- Eigen is great for **small matrices** and **expression templates**
- For **large dense linear algebra**, vendor-tuned BLAS/LAPACK (used by base R) is often faster
- The covariance computation + eigendecomposition is expensive and dominates runtime

## Solutions (In Priority Order)

### Solution 1: Use RcppArmadillo with ARPACK (RECOMMENDED)
Replace Eigen with Armadillo for whitening:
- Armadillo has `eigs_sym()` for computing only top-k eigenpairs
- Uses ARPACK under the hood (industry standard for sparse/truncated eigenproblems)
- Expected speedup: **5-10x faster whitening**

```cpp
// New file: inst/include/RcppICA/WhiteningArma.hpp
#include <RcppArmadillo.h>

// Compute only top k eigenvectors of covariance
arma::vec eigval;
arma::mat eigvec;
arma::eigs_sym(eigval, eigvec, C, k, "la");  // "la" = largest algebraic
```

### Solution 2: Implement Randomized PCA
For cases where n >> k, use randomized algorithms:
- Halko, Martinsson, Tropp (2011)
- Implemented in RcppML as `nmf_pca()`
- Expected speedup: **3-5x faster**

```cpp
// inst/include/RcppICA/RandomizedPCA.hpp
// Use power iteration for approximate top-k eigenvectors
```

### Solution 3: Call R's LAPACK Directly
Use R's built-in LAPACK via Rcpp:
```cpp
// Call R's dsyevr_ for symmetric eigendecomposition with range
// This is what base fastICA does
```

### Solution 4: Pre-compute Less
- Avoid forming explicit covariance matrix
- Use SVD of X_centered directly
- Eigen's BDCSVD, but with better flags

### Solution 5: Use OpenBLAS/MKL
Ensure Eigen is linked against optimized BLAS:
- Many speedups come from BLAS choice, not algorithm
- Eigen can use external BLAS via `EIGEN_USE_BLAS`
- But this requires user to have optimized BLAS installed

## Recommended Implementation Path

**Phase 1: Quick Fix (30 minutes)**
- [x] Changed default from SVD to Eigen (30% faster, already done)
- Update documentation to remove false speedup claims
- Add honest performance notes

**Phase 2: Real Fix (2-3 hours)**
- Add RcppArmadillo dependency
- Implement truncated eigendecomposition using `eigs_sym()`
- Expected result: 2-3x faster than base fastICA

**Phase 3: Optimization (Optional, 4-6 hours)**
- Implement randomized SVD for very large matrices
- Add heuristic to choose best method based on (n, m, k)
- Expected result: 5-10x faster than base fastICA

## Updated Performance Claims

### Current (After eigen default change)
- **Small data (n < 500)**: Similar speed to base fastICA (0.75x - 1.0x)
- **Medium data (n ≈ 1000)**: ~2x slower than base fastICA
- **Large data (n > 2000)**: ~3x slower than base fastICA

### After RcppArmadillo Fix (Projected)
- **Small data**: 1-2x faster than base fastICA
- **Medium data**: 2-3x faster than base fastICA
- **Large data**: 3-5x faster than base fastICA

### After Full Optimization (Projected)
- **Small data**: 2-3x faster
- **Medium data**: 5-10x faster
- **Large data**: 10-20x faster (with randomized methods)

## Code Changes Needed

### Immediate
1. Update vignettes to remove "30-150x speedup" claims
2. Add honest performance notes
3. Document that performance is currently competitive, not superior

### Short Term (RcppArmadillo Integration)
1. Add `LinkingTo: RcppArmadillo` to DESCRIPTION
2. Create new `inst/include/RcppICA/WhiteningArma.hpp`
3. Add conditional compilation: use Armadillo if available, else Eigen
4. Benchmark and verify speedup

### Medium Term
1. Implement randomized SVD
2. Add automatic method selection heuristic
3. Profile and optimize tight loops

## Why This Happened

Looking at the development history:
1. Package was ~85% complete when we started
2. Core algorithms were implemented WITHOUT benchmarking
3. Assumptions about Eigen performance were not validated
4. Focus was on correctness and tests, not performance
5. The "30-150x" claim appears to be aspirational, not measured

This is a good reminder to **always benchmark early** when performance is a key selling point.

## Next Steps

The user asked: "why would the user side be longer for us?"

**Answer**: Because our whitening implementation is inefficient. We compute full eigendecompositions when we only need partial ones. Base fastICA uses LAPACK more efficiently.

**Recommendation**: Implement Solution 1 (RcppArmadillo with ARPACK) before considering this package ready for release. Without it, we can't justify using RcppICA over base fastICA.
