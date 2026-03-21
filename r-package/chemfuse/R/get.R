#' Retrieve a single compound with optional enrichment
#'
#' Fetches detailed information for a compound by name or identifier and
#' optionally enriches the result with ADMET predictions, drug-likeness
#' filters, and binding data.
#'
#' @param query Character string. Compound name, CID, SMILES, or InChIKey.
#' @param admet Logical. Whether to include ADMET predictions. Defaults to
#'   \code{FALSE}.
#' @param druglikeness Logical. Whether to include drug-likeness filter
#'   results (Lipinski, Veber, PAINS). Defaults to \code{FALSE}.
#' @param binding Logical. Whether to retrieve binding data from BindingDB.
#'   Defaults to \code{FALSE}.
#' @param source Character string. Primary source database. Defaults to
#'   \code{"pubchem"}.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{compound}{Tibble with compound properties.}
#'     \item{admet}{Tibble with ADMET predictions, or \code{NULL}.}
#'     \item{druglikeness}{List of drug-likeness filter results, or \code{NULL}.}
#'     \item{binding}{Tibble with binding measurements, or \code{NULL}.}
#'   }
#'
#' @examples
#' \dontrun{
#'   # Basic retrieval
#'   result <- cf_get("aspirin")
#'
#'   # With ADMET predictions
#'   result <- cf_get("aspirin", admet = TRUE)
#'
#'   # Full enrichment
#'   result <- cf_get("aspirin", admet = TRUE, druglikeness = TRUE, binding = TRUE)
#' }
#'
#' @export
cf_get <- function(
    query,
    admet = FALSE,
    druglikeness = FALSE,
    binding = FALSE,
    source = "pubchem"
) {
  if (!is.character(query) || length(query) != 1L || nchar(query) == 0L) {
    stop("'query' must be a non-empty character string.")
  }

  # Search for the compound first
  collection <- cf$search(query, sources = list(source), limit = 1L)
  compounds <- collection$compounds

  if (length(compounds) == 0L) {
    return(list(
      compound    = tibble::tibble(),
      admet       = NULL,
      druglikeness = NULL,
      binding     = NULL
    ))
  }

  compound <- compounds[[1]]

  # Build the result list
  result <- list(
    compound = .compound_to_tibble(compound),
    admet = NULL,
    druglikeness = NULL,
    binding = NULL
  )

  # Optional ADMET enrichment
  if (isTRUE(admet) && !is.null(compound$smiles)) {
    tryCatch({
      pred <- cf$compute$admet$predict_admet(list(as.character(compound$smiles)))
      result$admet <- tibble::as_tibble(
        reticulate::py_to_r(pred) |> as.data.frame()
      )
    }, error = function(e) {
      message("ADMET prediction failed: ", conditionMessage(e))
    })
  }

  # Optional drug-likeness
  if (isTRUE(druglikeness) && !is.null(compound$smiles)) {
    tryCatch({
      from_compute <- reticulate::import("chemfuse.compute.druglikeness")
      filters <- from_compute$compute_all_filters(as.character(compound$smiles))
      result$druglikeness <- reticulate::py_to_r(filters)
    }, error = function(e) {
      message("Drug-likeness computation failed: ", conditionMessage(e))
    })
  }

  result
}


# Internal helper: convert a Python Compound to a single-row tibble
.compound_to_tibble <- function(compound) {
  props <- compound$properties

  tibble::tibble(
    cid             = if (is.null(compound$cid)) NA_character_ else as.character(compound$cid),
    name            = if (is.null(compound$name)) NA_character_ else as.character(compound$name),
    smiles          = if (is.null(compound$smiles)) NA_character_ else as.character(compound$smiles),
    inchikey        = if (is.null(compound$inchikey)) NA_character_ else as.character(compound$inchikey),
    inchi           = if (is.null(compound$inchi)) NA_character_ else as.character(compound$inchi),
    formula         = if (is.null(compound$formula)) NA_character_ else as.character(compound$formula),
    molecular_weight = if (is.null(props) || is.null(props$molecular_weight))
                          NA_real_ else as.numeric(props$molecular_weight),
    xlogp           = if (is.null(props) || is.null(props$xlogp))
                          NA_real_ else as.numeric(props$xlogp),
    hbd             = if (is.null(props) || is.null(props$hbd))
                          NA_integer_ else as.integer(props$hbd),
    hba             = if (is.null(props) || is.null(props$hba))
                          NA_integer_ else as.integer(props$hba),
    tpsa            = if (is.null(props) || is.null(props$tpsa))
                          NA_real_ else as.numeric(props$tpsa)
  )
}
