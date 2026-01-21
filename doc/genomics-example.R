## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)
library(RcppICA)
set.seed(42)

## ----simulate-data------------------------------------------------------------
# Simulate expression data similar to single-cell RNA-seq
n_cells <- 2000
n_genes <- 500

# Create independent gene programs
# Program 1: Cell cycle genes (periodic pattern)
program1 <- abs(sin(seq(0, 4*pi, length.out = n_cells)))

# Program 2: Immune response genes (binary activation)
program2 <- c(rep(1, 500), rep(0, 1000), rep(1, 500))

# Program 3: Metabolic genes (gradual increase)
program3 <- seq(0, 1, length.out = n_cells)^2

# Program 4: Differentiation trajectory (sigmoid)
program4 <- 1 / (1 + exp(-10 * (seq(-1, 1, length.out = n_cells))))

# Program 5: Housekeeping genes (constant with noise)
program5 <- rep(0.5, n_cells) + rnorm(n_cells, 0, 0.1)

# Combine into source matrix (cells × programs)
S_true <- cbind(program1, program2, program3, program4, program5)

# Create gene loading matrix (genes × programs)
# Each program affects different sets of genes
A_true <- matrix(0, n_genes, 5)

# Program 1 affects genes 1-100 (cell cycle)
A_true[1:100, 1] <- abs(rnorm(100, 1, 0.3))

# Program 2 affects genes 101-200 (immune)
A_true[101:200, 2] <- abs(rnorm(100, 1, 0.3))

# Program 3 affects genes 201-300 (metabolic)
A_true[201:300, 3] <- abs(rnorm(100, 1, 0.3))

# Program 4 affects genes 301-400 (differentiation)
A_true[301:400, 4] <- abs(rnorm(100, 1, 0.3))

# Program 5 affects genes 401-500 (housekeeping)
A_true[401:500, 5] <- abs(rnorm(100, 1, 0.3))

# Generate expression matrix: X = S * A^T
# Transpose to get genes × cells format (standard for genomics)
X_expression <- t(S_true %*% t(A_true))

# Add realistic noise
X_expression <- X_expression + abs(rnorm(length(X_expression), 0, 0.2))

# Log-transform (simulating log-normalized counts)
X_expression <- log1p(X_expression)

dim(X_expression)  # genes × cells

## ----run-ica------------------------------------------------------------------
# Transpose to cells × genes (ICA convention: observations × features)
X_cells <- t(X_expression)

# Run ICA with 5 components (matching our simulated programs)
cat("Running ICA on", nrow(X_cells), "cells and", ncol(X_cells), "genes...\n")

system.time({
  ica_result <- fastICA(X_cells,
                        n.comp = 5,
                        alg.typ = "parallel",  # Fast for multiple components
                        fun = "logcosh",       # Robust default
                        maxit = 300)
})

cat("\nConverged:", ica_result@misc$converged)
cat("\nIterations:", ica_result@misc$iterations)

## ----extract-results----------------------------------------------------------
# Independent components (ICs) - cell loadings (cells × components)
# Each row is a cell, each column is an IC
cell_loadings <- t(ica_result@S)
dim(cell_loadings)

# Mixing matrix - gene loadings (genes × components)
# Each row is a gene, each column shows which IC the gene contributes to
gene_loadings <- ica_result@A
dim(gene_loadings)

# Check recovery of original programs
# (Note: ICA may permute and flip signs)
recovery_cors <- abs(cor(cell_loadings, S_true))
cat("\nRecovery of true programs (best matches):\n")
print(round(apply(recovery_cors, 2, max), 3))

## ----interpret-programs-------------------------------------------------------
# For each IC, find top contributing genes
find_top_genes <- function(gene_loadings, ic_num, n_top = 20) {
  # Get absolute loadings for this IC
  loadings <- abs(gene_loadings[, ic_num])

  # Find top genes
  top_idx <- order(loadings, decreasing = TRUE)[1:n_top]

  data.frame(
    gene_id = top_idx,
    loading = gene_loadings[top_idx, ic_num],
    abs_loading = loadings[top_idx]
  )
}

# Example: Top genes for IC 1
top_genes_ic1 <- find_top_genes(gene_loadings, 1)
cat("\nTop 10 genes for IC 1:\n")
print(head(top_genes_ic1, 10))

# This IC should be enriched in genes 1-100 if it recovered the cell cycle program
cat("\nGene range 1-100 (cell cycle genes):",
    sum(top_genes_ic1$gene_id[1:20] <= 100), "out of 20 top genes\n")

## ----visualize-programs-------------------------------------------------------
# Plot cell loadings for each IC
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))

for (i in 1:5) {
  plot(cell_loadings[, i],
       type = "l",
       main = paste("Independent Component", i),
       xlab = "Cell Index",
       ylab = "Loading",
       col = rainbow(5)[i],
       lwd = 2)

  # Add reference line
  abline(h = 0, lty = 2, col = "gray50")
}

## ----compare-programs---------------------------------------------------------
# Find which IC corresponds to which true program
match_programs <- function(cell_loadings, S_true) {
  cors <- abs(cor(cell_loadings, S_true))

  # For each true program, find best matching IC
  matches <- apply(cors, 2, which.max)
  match_cors <- apply(cors, 2, max)

  data.frame(
    true_program = paste0("Program", 1:ncol(S_true)),
    matched_ic = paste0("IC", matches),
    correlation = round(match_cors, 3)
  )
}

program_matches <- match_programs(cell_loadings, S_true)
print(program_matches)

# Plot matched programs side by side
par(mfrow = c(5, 2), mar = c(3, 4, 2, 1))

for (i in 1:5) {
  matched_ic <- program_matches$matched_ic[i]
  ic_num <- as.numeric(sub("IC", "", matched_ic))

  # True program
  plot(S_true[, i], type = "l", col = "blue", lwd = 2,
       main = paste("True", program_matches$true_program[i]),
       ylab = "Signal", xlab = "")

  # Recovered IC (may need sign flip)
  ic_signal <- cell_loadings[, ic_num]
  if (cor(ic_signal, S_true[, i]) < 0) ic_signal <- -ic_signal

  plot(ic_signal, type = "l", col = "red", lwd = 2,
       main = paste("Recovered", matched_ic),
       ylab = "Signal", xlab = "")
}

## ----real-workflow, eval=FALSE------------------------------------------------
# # Example with SingleCellExperiment (not run)
# library(SingleCellExperiment)
# library(scuttle)
# 
# # Load your data
# # sce <- ZeiselBrainData()  # Example dataset
# 
# # Preprocessing
# # sce <- sce[!grepl("^mt-", rownames(sce)), ]  # Remove mitochondrial
# # sce <- logNormCounts(sce)                    # Log-normalize
# 
# # Extract expression matrix (genes × cells)
# # expr_matrix <- logcounts(sce)
# 
# # Transpose for ICA (cells × genes)
# # X <- t(as.matrix(expr_matrix))
# 
# # Feature selection (optional but recommended)
# # Select highly variable genes to reduce dimensionality
# # hvgs <- getTopHVGs(sce, n = 2000)
# # X_hvg <- X[, hvgs]
# 
# # Run ICA
# # ica_result <- fastICA(X_hvg,
# #                       n.comp = 15,           # Choose based on your data
# #                       alg.typ = "parallel",
# #                       maxit = 300,
# #                       seed = 42)
# 
# # Add to SingleCellExperiment
# # reducedDim(sce, "ICA") <- t(ica_result@S)
# 
# # Visualize on UMAP (after computing UMAP)
# # library(scater)
# # sce <- runUMAP(sce, dimred = "ICA")
# # plotReducedDim(sce, dimred = "UMAP", colour_by = "level1class")
# 
# # Plot specific IC on UMAP
# # ica_loadings <- t(ica_result@S)
# # plotUMAP(sce, colour_by = ica_loadings[, 1],
# #          point_size = 0.5) +
# #   ggtitle("ICA Component 1")

## ----benchmark-comparison-----------------------------------------------------
# Generate benchmark data
n_benchmark <- 1000
p_benchmark <- 500
X_benchmark <- matrix(rnorm(n_benchmark * p_benchmark), n_benchmark, p_benchmark)

# Time RcppICA
cat("RcppICA timing:\n")
time_rcpp <- system.time({
  result_rcpp <- fastICA(X_benchmark, n.comp = 10, maxit = 200)
})
print(time_rcpp)

# Compare with base fastICA if available
if (requireNamespace("fastICA", quietly = TRUE)) {
  cat("\nBase fastICA timing:\n")
  time_base <- system.time({
    result_base <- fastICA::fastICA(X_benchmark, n.comp = 10,
                                     maxit = 200, method = "C")
  })
  print(time_base)

  speedup <- time_base["elapsed"] / time_rcpp["elapsed"]
  cat("\nSpeedup:", round(speedup, 1), "×\n")

  # Verify results are similar (accounting for possible permutation)
  # RcppICA: S is (n_comp × n_obs), fastICA: S is (n_obs × n_comp)
  S_rcpp <- t(result_rcpp@S)  # Now (n_obs × n_comp)
  S_base <- result_base$S      # Already (n_obs × n_comp)

  cors <- abs(cor(S_rcpp, S_base))
  cat("Result correlation:", round(mean(apply(cors, 1, max)), 3), "\n")
} else {
  cat("\nfastICA package not available for comparison\n")
  cat("Install with: install.packages('fastICA')\n")
}

## ----feature-selection, eval=FALSE--------------------------------------------
# # Don't use all genes - select informative features
# # Option 1: Highly variable genes (recommended)
# # Use scran, Seurat, or similar to identify HVGs
# 
# # Option 2: Filter by expression
# # Remove genes expressed in < 5% of cells
# # expr_freq <- rowMeans(expr_matrix > 0)
# # informative_genes <- expr_freq > 0.05
# 
# # Option 3: Remove low-variance genes
# # gene_vars <- apply(expr_matrix, 1, var)
# # top_var_genes <- order(gene_vars, decreasing = TRUE)[1:2000]

## ----choosing-k, eval=FALSE---------------------------------------------------
# # Start with 10-20 components for exploration
# # ica_explore <- fastICA(X, n.comp = 15)
# 
# # Increase if you need finer resolution
# # ica_detailed <- fastICA(X, n.comp = 30)
# 
# # Unlike PCA, there's no variance explained metric
# # Evaluate based on:
# # - Biological interpretability
# # - Recovery of known cell types
# # - Pathway enrichment of top genes per IC

## ----interpretation-tips, eval=FALSE------------------------------------------
# # For each IC, perform gene set enrichment
# # library(enrichR)
# # top_genes <- find_top_genes(ica_result@A, ic_num = 1, n_top = 100)
# # enrichment <- enrichr(top_genes$gene_id, databases = "GO_Biological_Process_2021")
# 
# # Check if ICs correspond to known cell types
# # cor_with_markers <- cor(cell_loadings, known_marker_expression)
# 
# # Look for batch effects or technical artifacts
# # May appear as ICs that correlate with batch, depth, etc.

## ----postprocessing, eval=FALSE-----------------------------------------------
# # Rotate components if needed (ICA order is arbitrary)
# # order_by_variance <- order(apply(cell_loadings, 2, var), decreasing = TRUE)
# # cell_loadings_sorted <- cell_loadings[, order_by_variance]
# 
# # Remove technical components
# # Keep only biologically relevant ICs after inspection
# # bio_ics <- c(1, 2, 4, 5, 7)  # Example: exclude IC 3, 6 (technical)
# # cell_loadings_clean <- cell_loadings[, bio_ics]

## ----complete-example---------------------------------------------------------
# Simulate more realistic data with cell types
n_cells_per_type <- 400
n_cell_types <- 3
total_cells <- n_cells_per_type * n_cell_types

# Create cell type-specific programs
cell_type1 <- c(rep(1, n_cells_per_type), rep(0, n_cells_per_type * 2))
cell_type2 <- c(rep(0, n_cells_per_type), rep(1, n_cells_per_type),
                rep(0, n_cells_per_type))
cell_type3 <- c(rep(0, n_cells_per_type * 2), rep(1, n_cells_per_type))

# Add smooth noise for realism
cell_type1 <- cell_type1 + rnorm(total_cells, 0, 0.1)
cell_type2 <- cell_type2 + rnorm(total_cells, 0, 0.1)
cell_type3 <- cell_type3 + rnorm(total_cells, 0, 0.1)

# Shared housekeeping program
housekeeping <- rep(1, total_cells) + rnorm(total_cells, 0, 0.05)

# Combine
S_celltypes <- cbind(cell_type1, cell_type2, cell_type3, housekeeping)

# Gene loadings
A_celltypes <- matrix(0, 400, 4)
A_celltypes[1:100, 1] <- abs(rnorm(100, 2, 0.5))    # Type 1 markers
A_celltypes[101:200, 2] <- abs(rnorm(100, 2, 0.5))  # Type 2 markers
A_celltypes[201:300, 3] <- abs(rnorm(100, 2, 0.5))  # Type 3 markers
A_celltypes[301:400, 4] <- abs(rnorm(100, 1, 0.3))  # Housekeeping

# Generate expression: (cells × programs) × (programs × genes)^T = (cells × genes)
# Then transpose to (genes × cells)
X_celltypes_raw <- S_celltypes %*% t(A_celltypes)  # cells × genes
X_celltypes_raw <- X_celltypes_raw + matrix(rnorm(total_cells * 400, 0, 0.3), total_cells, 400)
X_celltypes <- t(X_celltypes_raw)  # Transpose to genes × cells, then back to cells × genes for ICA

# Run ICA
ica_celltypes <- fastICA(X_celltypes, n.comp = 4, maxit = 300)

# Visualize cell type separation
cell_loadings_ct <- t(ica_celltypes@S)

par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
for (i in 1:4) {
  plot(cell_loadings_ct[, i],
       col = rep(c("red", "blue", "green"), each = n_cells_per_type),
       pch = 16, cex = 0.5,
       main = paste("IC", i),
       xlab = "Cell Index",
       ylab = "Loading")
  abline(v = c(n_cells_per_type, n_cells_per_type * 2), lty = 2)
}

## ----session-info-------------------------------------------------------------
sessionInfo()

