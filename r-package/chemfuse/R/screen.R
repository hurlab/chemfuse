#' Batch screen compounds from a data frame or CSV file
#'
#' Searches for multiple compounds in parallel across chemical databases.
#' Accepts either a data frame with a \code{smiles} or \code{name} column,
#' or a path to a CSV file with those columns.
#'
#' @param data A data frame with at least one of \code{smiles} or \code{name}
#'   columns, or a character string giving the path to a CSV file.
#' @param sources Character vector. Database sources to query. Defaults to
#'   \code{c("pubchem")}.
#' @param admet Logical. Whether to compute ADMET predictions for each
#'   compound. Defaults to \code{FALSE}.
#' @param druglikeness Logical. Whether to compute drug-likeness filters.
#'   Defaults to \code{FALSE}.
#' @param concurrency Integer. Maximum concurrent requests. Defaults to 5.
#'
#' @return A tibble with one row per input compound containing the input
#'   identifiers plus retrieved database fields and optionally predicted
#'   ADMET properties.
#'
#' @examples
#' \dontrun{
#'   # Screen a data frame of SMILES
#'   df <- data.frame(smiles = c("CCO", "CCC", "CC(=O)O"))
#'   results <- cf_screen(df)
#'
#'   # Screen from a CSV file
#'   results <- cf_screen("compounds.csv", sources = c("pubchem", "chembl"))
#'
#'   # Include ADMET predictions
#'   results <- cf_screen(df, admet = TRUE)
#' }
#'
#' @export
cf_screen <- function(
    data,
    sources = c("pubchem"),
    admet = FALSE,
    druglikeness = FALSE,
    concurrency = 5L
) {
  # Resolve input to a data frame
  if (is.character(data) && length(data) == 1L) {
    if (!file.exists(data)) {
      stop("File not found: ", data)
    }
    data <- utils::read.csv(data, stringsAsFactors = FALSE)
  }

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame or a path to a CSV file.")
  }

  has_smiles <- "smiles" %in% names(data)
  has_name   <- "name"   %in% names(data)

  if (!has_smiles && !has_name) {
    stop("'data' must contain a 'smiles' or 'name' column.")
  }

  # Use the batch_search function from the Python package
  batch_fn <- reticulate::import("chemfuse.core.batch")

  # Write data to a temporary CSV for the Python batch function
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  utils::write.csv(data, tmp, row.names = FALSE)

  result <- cf$batch_search(
    tmp,
    sources = as.list(sources),
    concurrency = as.integer(concurrency)
  )

  collection <- result[[1]]
  .collection_to_tibble_extended(collection, data, admet, druglikeness)
}


# Internal helper: extend collection tibble with optional computed properties
.collection_to_tibble_extended <- function(collection, input_data, admet, druglikeness) {
  base_tbl <- .collection_to_tibble(collection)

  if (isTRUE(admet) && "smiles" %in% names(base_tbl) && nrow(base_tbl) > 0L) {
    smiles_list <- base_tbl$smiles[!is.na(base_tbl$smiles)]
    if (length(smiles_list) > 0L) {
      tryCatch({
        pred_mod <- reticulate::import("chemfuse.compute.admet")
        preds <- pred_mod$predict_admet(as.list(smiles_list))
        pred_df <- tibble::as_tibble(
          reticulate::py_to_r(preds) |> as.data.frame()
        )
        base_tbl <- dplyr::bind_cols(base_tbl, pred_df)
      }, error = function(e) {
        message("ADMET batch prediction failed: ", conditionMessage(e))
      })
    }
  }

  base_tbl
}
