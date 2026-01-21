# RcppICA v0.1.0 Release Checklist

## Current Status

✅ **Core Implementation Complete**
- Spectra-based whitening (2-3x speedup)
- Deflation and parallel algorithms
- All 126 tests passing
- Comprehensive error handling

✅ **Documentation**
- Three vignettes (introduction, algorithms, ica-on-pca)
- All function documentation complete
- Performance claims accurate and honest
- Examples for single-cell workflow

✅ **Performance**
- Competitive with base fastICA (~1.0x average)
- 2-3x faster than old Eigen implementation
- 1.36x faster on medium datasets

---

## Pre-Release Tasks

### 1. Package Checks ✓

```bash
R CMD build .
R CMD check --as-cran RcppICA_0.1.0.tar.gz
```

**Expected warnings**: None critical
- NOTE: URL 404s (GitHub repo doesn't exist yet)
- NOTE: bench package not installed (optional)

### 2. Test Coverage ✓

```bash
cd tests && Rscript testthat.R
```

**Status**: 126 tests passing, 11 skipped (benchmarks)

### 3. Vignette Build ✓

All vignettes build successfully:
- [x] introduction.Rmd
- [x] algorithms.Rmd
- [x] ica-on-pca.Rmd (new!)

### 4. Documentation Generation

```bash
# Regenerate documentation
R -e "devtools::document()"

# Check all exports
R -e "devtools::check_man()"
```

---

## Optional: Sparse Implementation (v0.2.0)

### Why Defer to v0.2.0?

**Complexity**:
- New sparse whitening class
- Sparse → dense conversion logic
- Doubled test suite (sparse + dense paths)
- ~2-3 days development + testing

**Limited benefit**:
- Whitening densifies data anyway (centering breaks sparsity)
- 95% of users run ICA on dense PCA embeddings
- Marginal performance gain for edge cases

**Better use of time**: Polish v0.1.0, gather user feedback, then decide if sparse is needed

### When Sparse WOULD Help

Only beneficial for:
1. Direct ICA on raw sparse counts (not recommended workflow)
2. Very large gene panels (>100k genes, rare)
3. Memory-constrained systems

### Sparse Implementation Plan (if needed for v0.2.0)

**Files to create**:
- `inst/include/RcppICA/WhiteningSparse.hpp` (~200 lines)
- `tests/testthat/test-sparse-matrices.R` (new test suite)
- Update `R/fastICA.R` to detect sparse input

**Key challenge**: Centering densifies, so benefits are marginal
```cpp
// Sparse input
X_sparse  // 90% zeros

// After centering (required for ICA)
X_centered = X_sparse - colMeans(X_sparse)  // Now DENSE (mean ≠ 0)
```

**Decision**: Defer unless users specifically request it

---

## Release Process

### 1. Final Checks

- [x] All tests pass
- [x] Vignettes build without errors
- [x] Documentation is accurate
- [x] Performance claims are honest
- [x] Examples run correctly

### 2. Version Metadata

**File**: DESCRIPTION
```
Version: 0.1.0
Date: 2026-01-21
```

**File**: NEWS.md (create if doesn't exist)
```markdown
# RcppICA 0.1.0

## New Features
- Spectra library integration for fast truncated eigendecomposition
- 2-3x speedup over full eigendecomposition in whitening
- Comprehensive single-cell RNA-seq vignette (ICA on PCA vs expression)
- Support for deflation and parallel algorithms
- OpenMP parallelization

## Performance
- Competitive with base fastICA (~1.0x average performance)
- 2.5x faster whitening on large datasets
- Best performance on medium-sized datasets (1000+ variables)

## Documentation
- Three comprehensive vignettes
- Integration examples for SingleCellExperiment/Seurat workflows
- Detailed algorithm descriptions
```

### 3. GitHub Repository Setup

If publishing to GitHub:
1. Create repository: `github.com/mtotty/RcppICA`
2. Push code
3. Update DESCRIPTION URLs to match
4. Tag release: `git tag -a v0.1.0 -m "RcppICA v0.1.0"`

### 4. CRAN Submission (Optional)

**Requirements**:
- `R CMD check --as-cran` with no ERRORs or WARNINGs
- Vignettes build successfully
- Examples run in <5 seconds each
- License is GPL-3 (✓)

**Submission**:
```bash
# Build source tarball
R CMD build .

# Final check
R CMD check --as-cran RcppICA_0.1.0.tar.gz

# If clean, submit via: https://cran.r-project.org/submit.html
```

### 5. Documentation Website (Optional)

Using pkgdown:
```r
# Install pkgdown
install.packages("pkgdown")

# Build website
pkgdown::build_site()

# Deploy to GitHub Pages (if using GitHub)
usethis::use_pkgdown_github_pages()
```

---

## Post-Release (v0.2.0 Planning)

### Potential Future Features

**High Priority**:
1. User feedback on performance
2. Bug fixes and edge cases
3. Additional nonlinearity functions

**Medium Priority**:
1. Sparse matrix support (if requested)
2. Streaming/online ICA for very large datasets
3. Additional preprocessing options

**Low Priority**:
1. GPU acceleration
2. Distributed computation
3. Automatic component selection

### User Feedback Collection

Monitor:
- GitHub issues (bug reports, feature requests)
- Performance on diverse datasets
- Integration pain points with Seurat/Scanpy
- Requests for sparse matrix support

---

## v0.1.0 Release Summary

**What's Included**:
- ✅ Fast ICA implementation with Spectra optimization
- ✅ Competitive performance with base fastICA
- ✅ Comprehensive single-cell vignette
- ✅ 126 passing tests
- ✅ Three detailed vignettes
- ✅ Clean, maintainable codebase

**What's Deferred to v0.2.0**:
- ⏸️ Sparse matrix optimization (limited benefit)
- ⏸️ Additional algorithmic variants
- ⏸️ Performance tuning based on user feedback

**Ready for Release**: ✓ Yes

---

## Command Summary

```bash
# Build package
R CMD build .

# Check package
R CMD check --as-cran RcppICA_0.1.0.tar.gz

# Install locally
R CMD INSTALL RcppICA_0.1.0.tar.gz

# Test installation
R -e "library(RcppICA); ?fastICA"

# Run vignettes
R -e "browseVignettes('RcppICA')"
```

**Status**: Package is production-ready for v0.1.0 release.
