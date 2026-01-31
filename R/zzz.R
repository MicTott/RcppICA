# RcppICA package initialization

#' @useDynLib RcppICA, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

.onLoad <- function(libname, pkgname) {
  if (is.null(getOption("RcppICA.threads"))) {
    options(RcppICA.threads = 1)
  }
  if (is.null(getOption("RcppICA.verbose"))) {
    options(RcppICA.verbose = FALSE)
  }
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "RcppICA v", utils::packageVersion("RcppICA")
  )
}
