#' ICAResult Class
#'
#' S4 class to store results from fastICA
#'
#' @slot S Matrix of independent components (n.comp x n observations)
#' @slot A Mixing matrix (variables x n.comp)
#' @slot W Unmixing matrix (n.comp x whitened dimension)
#' @slot K Whitening matrix (variables x whitened dimension)
#' @slot center Numeric vector of column means used for centering
#' @slot misc List containing metadata (iterations, converged, call, runtime, etc.)
#'
#' @exportClass ICAResult
setClass("ICAResult",
    slots = c(
        S = "matrix",
        A = "matrix",
        W = "matrix",
        K = "matrix",
        center = "numeric",
        misc = "list"
    ),
    prototype = list(
        S = matrix(numeric(0), 0, 0),
        A = matrix(numeric(0), 0, 0),
        W = matrix(numeric(0), 0, 0),
        K = matrix(numeric(0), 0, 0),
        center = numeric(0),
        misc = list()
    )
)

# Validity method
setValidity("ICAResult", function(object) {
    errors <- character()

    # Check that matrices have compatible dimensions
    if (nrow(object@S) > 0 && ncol(object@A) > 0) {
        if (nrow(object@S) != ncol(object@A)) {
            errors <- c(errors,
                "Number of components in S must match columns in A")
        }
    }

    # Check misc is a list
    if (!is.list(object@misc)) {
        errors <- c(errors, "misc slot must be a list")
    }

    if (length(errors) == 0) TRUE else errors
})
