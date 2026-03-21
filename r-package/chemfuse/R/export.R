#' Export results to CSV
#'
#' Writes a tibble or data frame of ChemFuse results to a CSV file.
#'
#' @param results A tibble or data frame of compound results as returned by
#'   \code{\link{cf_search}}, \code{\link{cf_screen}}, or similar functions.
#' @param path Character string. Output file path (e.g.,
#'   \code{"results.csv"}).
#' @param row.names Logical. Whether to include row names. Defaults to
#'   \code{FALSE}.
#'
#' @return Invisibly returns \code{path}.
#'
#' @examples
#' \dontrun{
#'   results <- cf_search("aspirin")
#'   cf_to_csv(results, "aspirin_results.csv")
#' }
#'
#' @export
cf_to_csv <- function(results, path, row.names = FALSE) {
  if (!is.data.frame(results)) {
    stop("'results' must be a data frame or tibble.")
  }
  if (!is.character(path) || length(path) != 1L || nchar(path) == 0L) {
    stop("'path' must be a non-empty character string.")
  }

  utils::write.csv(results, file = path, row.names = row.names)
  message("Results written to: ", path)
  invisible(path)
}


#' Export results to Excel
#'
#' Writes a tibble or data frame of ChemFuse results to an Excel (.xlsx)
#' file. Requires the \code{openxlsx} package to be installed.
#'
#' @param results A tibble or data frame of compound results as returned by
#'   \code{\link{cf_search}}, \code{\link{cf_screen}}, or similar functions.
#' @param path Character string. Output file path with \code{.xlsx}
#'   extension.
#' @param sheet Character string. Name of the worksheet. Defaults to
#'   \code{"compounds"}.
#'
#' @return Invisibly returns \code{path}.
#'
#' @examples
#' \dontrun{
#'   results <- cf_search("aspirin")
#'   cf_to_excel(results, "aspirin_results.xlsx")
#' }
#'
#' @export
cf_to_excel <- function(results, path, sheet = "compounds") {
  if (!is.data.frame(results)) {
    stop("'results' must be a data frame or tibble.")
  }
  if (!is.character(path) || length(path) != 1L || nchar(path) == 0L) {
    stop("'path' must be a non-empty character string.")
  }
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop(
      "Package 'openxlsx' is required for Excel export. ",
      "Install it with: install.packages('openxlsx')"
    )
  }

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheet)
  openxlsx::writeData(wb, sheet, results)
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)

  message("Results written to: ", path)
  invisible(path)
}
