# RcppICA package initialization
#
# This file contains package initialization functions that run when
# the package is loaded or attached.

#' @useDynLib RcppICA, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

.onLoad <- function(libname, pkgname) {
  # Set default package options if not already set
  if (is.null(getOption("RcppICA.threads"))) {
    options(RcppICA.threads = 0)  # 0 = auto-detect all available cores
  }

  if (is.null(getOption("RcppICA.verbose"))) {
    options(RcppICA.verbose = FALSE)
  }
}

.onAttach <- function(libname, pkgname) {
  # Get thread configuration
  threads <- getOption("RcppICA.threads", default = 0)

  # Determine actual thread count for display
  if (threads == 0) {
    # Auto-detect
    actual_threads <- parallel::detectCores(logical = FALSE)
    if (is.na(actual_threads)) {
      actual_threads <- 1
    }
    thread_msg <- sprintf("Using %d threads (auto-detected)", actual_threads)
  } else {
    thread_msg <- sprintf("Using %d thread(s)", threads)
  }

  # Package startup message
  packageStartupMessage(
    "RcppICA v", utils::packageVersion("RcppICA"), "\n",
    "Ultra-fast Independent Component Analysis using Rcpp and Eigen\n",
    thread_msg, "\n",
    "Set options: options(RcppICA.threads = N, RcppICA.verbose = TRUE/FALSE)"
  )
}
