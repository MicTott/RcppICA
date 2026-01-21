# RcppICA

<!-- badges: start -->
[![R-CMD-check](https://github.com/mtotty/RcppICA/workflows/R-CMD-check/badge.svg)](https://github.com/mtotty/RcppICA/actions)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Project Status: Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
<!-- badges: end -->

Fast Independent Component Analysis (ICA) using C++ with Eigen and Spectra libraries.

## Overview

**RcppICA** provides a high-performance implementation of the FastICA algorithm using optimized C++ code. It features:

- 🚀 **Fast truncated eigendecomposition** using the Spectra library (2-3x speedup)
- ⚡ **Competitive performance** with base fastICA (~1.0x average, faster on medium datasets)
- 🧬 **Single-cell RNA-seq ready** with comprehensive vignettes and workflows
- 🔧 **Flexible algorithms**: Parallel (symmetric) and deflation modes
- 🧵 **OpenMP parallelization** for multi-core systems
- 📊 **S4 classes** for clean object-oriented interface

## Installation

### From GitHub

```r
# Install devtools if needed
install.packages("devtools")

# Install RcppICA
devtools::install_github("mtotty/RcppICA")
```

### From Source

```bash
# Clone repository
git clone https://github.com/mtotty/RcppICA.git
cd RcppICA

# Build and install
R CMD build .
R CMD INSTALL RcppICA_0.1.0.tar.gz
```

### System Requirements

- R (>= 4.0.0)
- C++17 compiler
- GNU make

## Quick Start

### Basic Usage

```r
library(RcppICA)

# Generate mixed signals (cocktail party problem)
set.seed(42)
n <- 1000

# Two independent sources
S <- cbind(
  sin(seq(0, 8*pi, length.out = n)),     # Sinusoid
  ((1:n %% 50) - 25) / 25                 # Sawtooth
)

# Random mixing
A <- matrix(c(0.5, 0.8, 0.6, 0.3), 2, 2)
X <- S %*% A

# Recover independent components
result <- fastICA(X, n.comp = 2)

# Check convergence
print(result)

# Compare recovered vs original (allowing for sign/order)
cor(t(result@S), S)
```

### Single-Cell RNA-seq Workflow

**Recommended approach**: Run ICA on PCA components

```r
library(RcppICA)
library(SingleCellExperiment)
library(scater)

# Standard preprocessing
sce <- logNormCounts(sce)
sce <- runPCA(sce, ncomponents = 50)

# Run ICA on top 30 PCs
ica_result <- fastICA(
    reducedDim(sce, "PCA")[, 1:30],  # cells × 30 PCs
    n.comp = 15,                      # Extract 15 ICs
    whiten.method = "spectra"         # Fast (default)
)

# Store in SCE object
reducedDim(sce, "ICA") <- t(ica_result@S)

# Calculate gene weights by projection
pca_loadings <- attr(reducedDim(sce, "PCA"), "rotation")[, 1:30]
gene_weights <- pca_loadings %*% ica_result@A
# Result: genes × 15 ICs
```

See the **"ICA on PCA vs Expression"** vignette for detailed comparison and best practices.

## Performance

Benchmark results on representative datasets:

| Dataset | Base fastICA | RcppICA (Spectra) | Speedup |
|---------|--------------|-------------------|---------|
| 500×100, k=5 | 0.008s | 0.009s | 0.89x |
| **1000×500, k=5** | 0.045s | **0.033s** | **1.36x ✓** |
| 5000×1000, k=10 | 0.363s | 0.501s | 0.72x |

**Average**: ~1.0x competitive performance across dataset sizes

**Whitening speedup**: 2.5x faster than full eigendecomposition (internal)

### Why Competitive, Not Faster?

Base fastICA uses highly optimized LAPACK (Accelerate/MKL). RcppICA's advantage:
- **Algorithmic improvement**: Computes only top-k eigenvalues (not all m)
- **Best on medium datasets**: 1.36x faster where k << m advantage is strongest
- **Cleaner code**: Modern C++, maintainable, S4 classes
- **Room for optimization**: Future versions can tune further

## Features

### Whitening Methods

- `"spectra"` (default): Fast truncated eigendecomposition using Lanczos iteration
- `"eigen"`: Full eigendecomposition of covariance matrix
- `"svd"`: SVD-based (most numerically stable)

### Algorithms

- `"parallel"` (default): Symmetric orthogonalization, faster and more stable
- `"deflation"`: Sequential extraction, may be preferred for few components

### Nonlinearity Functions

- `"logcosh"` (default): Most robust, works for sub- and super-Gaussian
- `"exp"`: Good for super-Gaussian distributions
- `"cube"`: Kurtosis-based, simplest

## Documentation

### Vignettes

- **Introduction to RcppICA**: Basic usage and examples
- **Algorithm Details**: FastICA theory and implementation
- **ICA on PCA vs Expression**: Single-cell best practices

View vignettes:
```r
browseVignettes("RcppICA")
```

### Help

```r
?fastICA           # Main function documentation
?ICAResult-class   # S4 result class
```

## Single-Cell RNA-seq: Best Practices

### ✓ Recommended: ICA on PCA Components

```r
# 1. Run PCA first (standard scRNA-seq workflow)
sce <- runPCA(sce, ncomponents = 50)

# 2. Run ICA on PCA embeddings
ica_result <- fastICA(reducedDim(sce, "PCA")[, 1:30], n.comp = 15)

# 3. Get gene weights via projection
gene_weights <- pca_loadings %*% ica_result@A
```

**Why?**
- 5-10x faster (30 PCs vs 2000 genes)
- Filters technical noise
- Standard field practice
- Same biological insights

### ✗ Not Recommended: ICA on Full Expression

```r
# Slower, higher memory, includes noise
ica_result <- fastICA(t(logcounts(sce)), n.comp = 15)
```

**Only use for**: Small gene panels (<500 genes)

See the **ica-on-pca** vignette for comprehensive comparison.

## Integration

### Seurat

```r
library(Seurat)
library(RcppICA)

# Run ICA on PCA
pca_embeddings <- Embeddings(seurat_obj, "pca")[, 1:30]
ica_result <- fastICA(pca_embeddings, n.comp = 15)

# Calculate gene weights
pca_loadings <- Loadings(seurat_obj, "pca")[, 1:30]
gene_weights <- pca_loadings %*% ica_result@A

# Store in Seurat
seurat_obj[["ica"]] <- CreateDimReducObject(
    embeddings = t(ica_result@S),
    loadings = gene_weights,
    key = "IC_"
)
```

### SingleCellExperiment

```r
library(SingleCellExperiment)

# Run ICA on PCA
ica_result <- fastICA(reducedDim(sce, "PCA")[, 1:30], n.comp = 15)

# Store results
reducedDim(sce, "ICA") <- t(ica_result@S)
```

## Examples

### Identifying Cell-Type Programs

```r
# Get top genes for each IC
get_top_genes <- function(gene_weights, ic, n = 20) {
    weights <- gene_weights[, ic]
    sort(abs(weights), decreasing = TRUE)[1:n]
}

top_genes_ic1 <- get_top_genes(gene_weights, ic = 1)
print(names(top_genes_ic1))
```

### Pathway Enrichment

```r
library(clusterProfiler)

# Top genes for IC 1
ic1_genes <- names(get_top_genes(gene_weights, 1, n = 200))

# GO enrichment
ego <- enrichGO(
    gene = ic1_genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP"
)
```

## Development

### Building from Source

```bash
# Clone repository
git clone https://github.com/mtotty/RcppICA.git
cd RcppICA

# Install dependencies
R -e "install.packages(c('Rcpp', 'RcppEigen', 'Matrix', 'testthat'))"

# Build
R CMD build .

# Check
R CMD check --as-cran RcppICA_0.1.0.tar.gz

# Install
R CMD INSTALL RcppICA_0.1.0.tar.gz
```

### Running Tests

```r
# Run all tests
devtools::test()

# Run specific test file
testthat::test_file("tests/testthat/test-fastICA.R")
```

### Benchmarking

```r
# Run included benchmarks
source("benchmark_spectra.R")
source("profile_spectra.R")
```

## Technical Details

### Implementation

- **C++ core**: Eigen 3.4+ for linear algebra
- **Spectra library**: v1.0.1 for truncated eigendecomposition (header-only)
- **OpenMP**: Parallel algorithm uses multi-threading
- **S4 classes**: Clean R interface with `ICAResult` class

### Algorithms

FastICA maximizes non-Gaussianity to find independent components:

1. **Whitening**: Decorrelate data and scale to unit variance
   - Spectra: Computes only top-k eigenvalues (fast)
   - Eigen: Full eigendecomposition (slower)
   - SVD: Most numerically stable

2. **Fixed-point iteration**: Optimize non-Gaussianity measure
   - Parallel: Symmetric orthogonalization (default)
   - Deflation: Gram-Schmidt orthogonalization

3. **Convergence**: Stops when unmixing matrix stabilizes

### Complexity

- **Spectra whitening**: O(n·m² + k·m²·iter) where iter ≈ 20-50
- **Full eigendecomposition**: O(n·m² + m³)
- **Speedup**: ~m/k theoretical (e.g., 500/5 = 100x), ~2.5x practical

## Citation

If you use RcppICA in published research, please cite:

```
Totty M (2026). RcppICA: Fast Independent Component Analysis Using Rcpp and Eigen.
R package version 0.1.0. https://github.com/mtotty/RcppICA
```

**FastICA algorithm**:
```
Hyvärinen A, Oja E (2000). Independent Component Analysis: Algorithms and
Applications. Neural Networks, 13(4-5):411-430.
```

## License

GPL-3

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure `R CMD check` passes
5. Submit a pull request

## Issues

Report bugs and request features at: https://github.com/mtotty/RcppICA/issues

## Acknowledgments

- **Spectra library**: Yixuan Qiu (https://spectralib.org/)
- **Eigen library**: Eigen development team
- **FastICA algorithm**: Ari Hyvärinen and Erkki Oja

---

**RcppICA** is actively maintained. For questions, please open an issue on GitHub.
