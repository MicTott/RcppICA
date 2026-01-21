#' @include AllClasses.R
NULL

#' Show Method for ICAResult
#'
#' @param object An ICAResult object
#' @export
setMethod("show", "ICAResult", function(object) {
    cat("ICA Result\n")
    cat("-----------\n")
    cat("Components:", nrow(object@S), "\n")
    cat("Observations:", ncol(object@S), "\n")
    cat("Variables:", nrow(object@A), "\n")
    if (!is.null(object@misc$iterations)) {
        cat("Iterations:", object@misc$iterations, "\n")
    }
    if (!is.null(object@misc$converged)) {
        cat("Converged:", object@misc$converged, "\n")
    }
})

#' Summary Method for ICAResult
#'
#' @param object An ICAResult object
#' @return Invisibly returns the object
#' @export
setMethod("summary", "ICAResult", function(object) {
    cat("Independent Component Analysis Results\n")
    cat("======================================\n\n")
    if (!is.null(object@misc$call)) {
        cat("Call:\n")
        print(object@misc$call)
    }
    cat("\nDimensions:\n")
    cat("  Independent components (S):",
        paste(dim(object@S), collapse = " x "), "\n")
    cat("  Mixing matrix (A):",
        paste(dim(object@A), collapse = " x "), "\n")
    cat("  Unmixing matrix (W):",
        paste(dim(object@W), collapse = " x "), "\n")
    cat("  Whitening matrix (K):",
        paste(dim(object@K), collapse = " x "), "\n")
    cat("\nConvergence:\n")
    if (!is.null(object@misc$iterations)) {
        cat("  Iterations:", object@misc$iterations, "\n")
    }
    if (!is.null(object@misc$converged)) {
        cat("  Converged:", object@misc$converged, "\n")
    }
    if (!is.null(object@misc$runtime)) {
        cat("  Runtime:", round(as.numeric(object@misc$runtime), 3), "secs\n")
    }

    # Component statistics
    cat("\nComponent Statistics:\n")
    cat("  Mean of |S|:", round(mean(abs(object@S)), 4), "\n")
    cat("  SD of S:", round(sd(as.vector(object@S)), 4), "\n")

    invisible(object)
})

#' Predict Method for ICAResult
#'
#' Project new data onto learned independent components
#'
#' @param object An ICAResult object
#' @param newdata Optional new data matrix (observations x variables).
#'   If missing, returns transposed S from original fit.
#' @return Matrix of independent component scores (observations x n.comp)
#'
#' @examples
#' \dontrun{
#' # Fit ICA
#' result <- fastICA(X_train, n.comp = 5)
#'
#' # Project new data
#' S_new <- predict(result, X_test)
#' }
#'
#' @export
setMethod("predict", "ICAResult", function(object, newdata) {
    if (missing(newdata)) {
        return(t(object@S))
    }

    # Validate newdata
    if (!is.matrix(newdata)) {
        newdata <- as.matrix(newdata)
    }

    if (ncol(newdata) != length(object@center)) {
        stop("newdata must have the same number of variables as original data (",
             length(object@center), " variables)")
    }

    # Center using stored means
    newdata_centered <- sweep(newdata, 2, object@center)

    # Project: S_new = (X_new - mean) * K * W^T
    S_new <- newdata_centered %*% object@K %*% t(object@W)

    return(S_new)
})

#' Plot Method for ICAResult
#'
#' Plot independent components or mixing matrix
#'
#' @param x An ICAResult object
#' @param y Ignored (for S4 compatibility)
#' @param type Type of plot: \code{"components"} (default) shows each IC
#'   as a time series, \code{"mixing"} shows heatmap of mixing matrix
#' @param components Which components to plot (default: all)
#' @param ... Additional arguments passed to plot functions
#'
#' @examples
#' \dontrun{
#' result <- fastICA(X, n.comp = 3)
#'
#' # Plot all components
#' plot(result)
#'
#' # Plot specific components
#' plot(result, components = 1:2)
#'
#' # Plot mixing matrix
#' plot(result, type = "mixing")
#' }
#'
#' @export
setMethod("plot", signature(x = "ICAResult", y = "missing"),
    function(x, y, type = c("components", "mixing"), components = NULL, ...) {
        type <- match.arg(type)

        if (type == "components") {
            n_comp <- nrow(x@S)

            if (is.null(components)) {
                components <- seq_len(n_comp)
            }

            components <- as.integer(components)
            components <- components[components >= 1 & components <= n_comp]

            if (length(components) == 0) {
                stop("No valid components to plot")
            }

            n_plots <- length(components)
            old_par <- par(mfrow = c(n_plots, 1), mar = c(2, 4, 2, 1))
            on.exit(par(old_par))

            for (i in seq_along(components)) {
                comp_idx <- components[i]
                plot(x@S[comp_idx, ], type = "l",
                     ylab = paste("IC", comp_idx),
                     main = if(i == 1) "Independent Components" else "",
                     xlab = if(i == n_plots) "Observation" else "",
                     ...)
            }
        } else {
            # Heatmap of mixing matrix
            old_par <- par(mar = c(5, 4, 4, 2) + 0.1)
            on.exit(par(old_par))

            image(t(x@A)[, nrow(x@A):1, drop = FALSE],
                  main = "Mixing Matrix (A)",
                  xlab = "Components",
                  ylab = "Variables",
                  axes = FALSE,
                  col = hcl.colors(100, "RdBu", rev = TRUE),
                  ...)

            # Add axes
            n_comp <- ncol(x@A)
            n_var <- nrow(x@A)
            axis(1, at = seq(0, 1, length.out = n_comp),
                 labels = seq_len(n_comp))
            axis(2, at = seq(0, 1, length.out = min(10, n_var)),
                 labels = round(seq(n_var, 1, length.out = min(10, n_var))))
        }

        invisible(x)
    }
)

#' Extract Independent Components
#'
#' @param object An ICAResult object
#' @return Matrix of independent components (n.comp x n observations)
#' @export
setGeneric("components", function(object) standardGeneric("components"))

#' @rdname components
#' @export
setMethod("components", "ICAResult", function(object) {
    object@S
})

#' Extract Mixing Matrix
#'
#' @param object An ICAResult object
#' @return Mixing matrix (variables x n.comp)
#' @export
setGeneric("mixing", function(object) standardGeneric("mixing"))

#' @rdname mixing
#' @export
setMethod("mixing", "ICAResult", function(object) {
    object@A
})
