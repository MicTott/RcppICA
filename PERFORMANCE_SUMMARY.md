# RcppICA Performance Investigation Summary

## Your Question

> "why would the user side be longer for us?"

## Answer

**RcppICA is currently 2-3x SLOWER than base fastICA**, not faster. The "user time" is longer because our whitening (preprocessing) step is inefficient.

## Root Cause

The bottleneck is in eigendecomposition-based whitening ([Whitening.hpp:90-94](inst/include/RcppICA/Whitening.hpp#L90-L94)):

```cpp
// We compute the FULL covariance matrix
Matrix C;
C.noalias() = X_centered.transpose() * X_centered;  // m × m matrix
C /= static_cast<Scalar>(n - 1);

// We compute ALL eigenvalues/eigenvectors
Eigen::SelfAdjointEigenSolver<Matrix> eig(C);  // All m eigenpairs
```

**The problem**: For a 1000-variable dataset, we:
1. Compute all 1000 eigenvalues/eigenvectors
2. Use only the top 5-10
3. Discard 99% of the computation

Base fastICA uses R's highly-optimized LAPACK more efficiently.

## Benchmark Results

| Dataset Size | RcppICA (eigen) | Base fast ICA (C) | Performance |
|--------------|-----------------|-------------------|-------------|
| 500×100, k=5 | 0.012s | 0.009s | **1.33x SLOWER** |
| 1000×500, k=5 | 0.105s | 0.047s | **2.23x SLOWER** |
| 5000×1000, k=10 | 1.267s | 0.426s | **2.97x SLOWER** |

The performance gap **increases with data size**, confirming algorithmic inefficiency.

## Time Breakdown

Whitening consumes 95-98% of total runtime:

```
Total RcppICA time (1000×500):  0.105s
  ├─ Whitening:                 0.103s (98%)
  └─ ICA iterations:            0.002s (2%)
```

**The ICA algorithm itself is very fast** - the problem is purely preprocessing.

## Why The Original Claim Was Wrong

The package plan claimed "30-150x speedup" based on:
- **Incorrect assumption**: That Eigen would be faster than R's LAPACK
- **No benchmarking**: Claims were aspirational, not measured
- **Wrong focus**: ICA iterations ARE fast, but they're not the bottleneck

**Reality**:
- Eigen is great for small matrices and expression templates
- For large dense linear algebra, vendor-tuned BLAS/LAPACK (used by R) is faster
- The covariance + eigendecomposition dominates runtime

## Attempted Fix: RcppArmadillo

I attempted to add RcppArmadillo for truncated eigendecomposition using ARPACK:
- Armadillo's `eigs_sym()` computes only top-k eigenpairs
- Expected 5-10x faster whitening

**Result**: FAILED due to header conflicts
- RcppArmadillo requires NOT including `Rcpp.h`
- RcppExports.cpp (auto-generated) includes both `Rcpp.h` and `RcppArmadillo.h`
- Error: "The file 'Rcpp.h' should not be included..."

**Solutions exist** (manual exports, separate compilation units) but are complex and time-consuming.

## Current State

**Package Status**:
- ✓ Builds and installs successfully
- ✓ All 165 tests pass
- ✓ Algorithms are correct
- ✗ Performance is 2-3x SLOWER than base fastICA

**Default Settings**:
- Whitening: `eigen` (30% faster than SVD, but still slow)
- Algorithm: `parallel` (good)
- Nonlinearity: `logcosh` (good)

## Recommendations

### Option 1: Document Honestly (IMMEDIATE)
**Time**: 30 minutes
**Impact**: Manage expectations

- Remove "30-150x speedup" claims from all documentation
- Add honest performance notes to vignettes
- State that RcppICA is "competitive" but not faster
- Focus on other benefits: S4 classes, modern interface, extensive testing

### Option 2: Implement Truncated Eigendecomposition (SHORT TERM)
**Time**: 4-6 hours
**Impact**: 2-3x faster than base fastICA

Approaches:
1. **RcppArmadillo** (best, but complex due to header conflicts)
2. **Call R's eigen() from C++** (simple but has R call overhead)
3. **Implement randomized SVD** (good for k << m)
4. **Use Spectra library** (C++ ARPACK wrapper, clean)

### Option 3: Accept Current Performance (PRAGMATIC)
**Time**: None
**Impact**: Package remains usable but not superior

Arguments for:
- Package is correct and well-tested
- Performance is "competitive" (~2x slower isn't terrible)
- Base fastICA is already highly optimized C code
- Focus on other value: maintainability, modern R interface, S4 classes

## My Recommendation

**For immediate release**: Go with Option 1 + Option 3
- Update all documentation to be honest about performance
- Position RcppICA as a "modern, well-tested alternative" not a "faster replacement"
- Highlight real benefits: S4 classes, extensive documentation, comprehensive tests

**For future v0.2.0**: Implement Option 2
- Add Spectra library for truncated eigendecomposition
- Achieve genuine 2-5x speedup over base fastICA
- Then you can legitimately claim performance benefits

## Files Created

1. **[PERFORMANCE_ANALYSIS.md](PERFORMANCE_ANALYSIS.md)** - Detailed technical analysis
2. **[fair_benchmark.R](fair_benchmark.R)** - Comprehensive benchmarking script
3. **[profile_performance.R](profile_performance.R)** - Profiling script
4. **[test_whitening_performance.R](test_whitening_performance.R)** - Whitening benchmarks

## Next Steps

If you want me to:
1. **Update documentation honestly** - I can remove false claims and add realistic performance notes
2. **Implement fast whitening** - I can add Spectra or another solution (4-6 hours)
3. **Do nothing** - Package works correctly as-is, just slower than claimed

Let me know which direction you'd like to go.
