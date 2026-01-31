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
    options(RcppICA.threads = 1)  # Default to 1 thread (safe for HPC)
  }

  if (is.null(getOption("RcppICA.verbose"))) {
    options(RcppICA.verbose = FALSE)
  }
}

.onAttach <- function(libname, pkgname) {
  # Get thread configuration
  threads <- getOption("RcppICA.threads", default = 1)

  # Build thread message

  thread_msg <- sprintf("Using %d thread(s)", threads)
  if (threads == 1) {
    thread_msg <- paste0(thread_msg, " (set n.threads or options(RcppICA.threads=N) to parallelize)")
  }

  # Package startup message
  packageStartupMessage(
    "RcppICA v", utils::packageVersion("RcppICA"), "\n",
    "Ultra-fast Independent Component Analysis using Rcpp and Eigen\n",
    thread_msg
  )
}
