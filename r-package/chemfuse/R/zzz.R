# Package-level environment for the Python module reference
cf <- NULL

#' @importFrom reticulate import
.onLoad <- function(libname, pkgname) {
  cf <<- reticulate::import("chemfuse", delay_load = TRUE)
}
