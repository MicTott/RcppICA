# RcppICA

<!-- badges: start -->
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

Fast Independent Component Analysis (ICA) for R, implemented in C++ with Eigen, Spectra, and OpenMP.

## Overview

**RcppICA** implements the FastICA algorithm with:

- **BLAS-accelerated linear algebra** via Eigen with vendor BLAS delegation (Accelerate/MKL/OpenBLAS)
- **Truncated eigendecomposition** using the Spectra library (computes only top-k eigenvalues)
- **Memory-efficient whitening** that avoids materializing intermediate n x m matrices
- **Sparse matrix support** for large single-cell datasets
- **OpenMP parallelization** for the ICA iteration
- **S4 classes** for a clean R interface

## Installation

```r
# Install from GitHub
devtools::install_github("MicTott/RcppICA")
```

### System Requirements

- R (>= 4.0.0)
- C++17 compiler

## Quick Start

```r
library(RcppICA)

# Generate mixed signals
set.seed(42)
n <- 1000
S <- cbind(
  sin(seq(0, 8*pi, length.out = n)),
  ((1:n %% 50) - 25) / 25
)
A <- matrix(c(0.5, 0.8, 0.6, 0.3), 2, 2)
X <- S %*% A

# Recover independent components
result <- fastICA(X, n.comp = 2)

# Check recovery (allowing for sign/order ambiguity)
cor(t(result@S), S)
```

### With Sparse Matrices

```r
library(Matrix)

# Pass sparse matrix directly - optimal path is selected automatically
X_sparse <- as(X, "dgCMatrix")
result <- fastICA(X_sparse, n.comp = 10)
```

## Performance

Benchmarks on an Apple M1 (10 cores), comparing RcppICA against the `fastICA` package:

| Dataset | RcppICA | fastICA | Speedup |
|---------|---------|---------|---------|
| 1k x 500, k=10 | 0.06s | 0.12s | 2x |
| 5k x 2k, k=10 | 0.7s | 2.4s | 3x |
| 10k x 5k, k=10 | 4.9s | 33s | 7x |
| 10k x 10k, k=10 | 17s | 230s | 14x |

The speedup increases with matrix size due to memory-efficient whitening and BLAS delegation.

## Features

### Whitening Methods

- `"spectra"` (default): Truncated eigendecomposition via Lanczos iteration
- `"eigen"`: Full eigendecomposition of covariance matrix
- `"svd"`: SVD-based (most numerically stable, uses covariance approach for large matrices)

### Algorithms

- `"parallel"` (default): Symmetric orthogonalization, extracts all components simultaneously
- `"deflation"`: Sequential extraction via Gram-Schmidt

### Nonlinearity Functions

- `"logcosh"` (default): Robust, works for sub- and super-Gaussian
- `"exp"`: Good for super-Gaussian (heavy-tailed) distributions
- `"cube"`: Kurtosis-based, simplest

## Single-Cell RNA-seq Usage

Run ICA on PCA embeddings (recommended):

```r
library(SingleCellExperiment)

# Standard preprocessing
sce <- logNormCounts(sce)
sce <- runPCA(sce, ncomponents = 50)


# Sparse path avoids full matrix materialization
ica_result <- fastICA(t(logcounts(sce)), n.comp = 20)
```

## Documentation

```r
browseVignettes("RcppICA")
?fastICA
?ICAResult-class
```

## License

GPL-3

## Acknowledgments

This package was heavily inspired by and modeled around the NMF implementation in [RcppML](https://github.com/zdebruine/RcppML).

## References

- Hyvarinen A, Oja E (2000). Independent Component Analysis: Algorithms and Applications. *Neural Networks*.
- DeBruine ZJ, Melber K, Bhatt DK, et al. (2021). Fast and robust non-negative matrix factorization for single-cell experiments. *bioRxiv*.
- Spectra library: https://spectralib.org/
