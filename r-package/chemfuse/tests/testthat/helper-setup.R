# Helper: check if Python chemfuse is available for integration tests
has_python_chemfuse <- function() {
  tryCatch({
    reticulate::import("chemfuse")
    TRUE
  }, error = function(e) {
    FALSE
  })
}

# Skip integration tests if Python chemfuse is not installed
skip_if_no_python <- function() {
  testthat::skip_if_not(
    has_python_chemfuse(),
    message = "Python 'chemfuse' package not available"
  )
}
