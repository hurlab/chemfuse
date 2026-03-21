#' Search for compounds across chemical databases
#'
#' Searches one or more chemical databases (PubChem, ChEMBL, UniChem, etc.)
#' and returns a tibble of matching compounds with merged, deduplicated results.
#'
#' @param query Character string. The search term (compound name, SMILES,
#'   CID, formula, or InChI).
#' @param sources Character vector. Database sources to query. Defaults to
#'   \code{c("pubchem")}. Valid values include \code{"pubchem"},
#'   \code{"chembl"}, \code{"unichem"}, \code{"bindingdb"},
#'   \code{"opentargets"}, \code{"surechembl"}.
#' @param query_type Character string. Type of query: \code{"name"},
#'   \code{"smiles"}, \code{"cid"}, \code{"formula"}, or \code{"inchi"}.
#'   Defaults to \code{"name"}.
#' @param limit Integer. Maximum results per source. Defaults to 100.
#'
#' @return A tibble with one row per compound and columns including
#'   \code{cid}, \code{name}, \code{smiles}, \code{inchikey}, \code{formula},
#'   \code{molecular_weight}, and \code{sources}.
#'
#' @examples
#' \dontrun{
#'   # Search by name
#'   results <- cf_search("aspirin")
#'
#'   # Search multiple sources
#'   results <- cf_search("aspirin", sources = c("pubchem", "chembl"))
#'
#'   # Search by SMILES
#'   results <- cf_search("CC(=O)Oc1ccccc1C(=O)O", query_type = "smiles")
#' }
#'
#' @export
cf_search <- function(
    query,
    sources = c("pubchem"),
    query_type = "name",
    limit = 100L
) {
  if (!is.character(query) || length(query) != 1L || nchar(query) == 0L) {
    stop("'query' must be a non-empty character string.")
  }
  if (!is.character(sources) || length(sources) == 0L) {
    stop("'sources' must be a non-empty character vector.")
  }

  result <- cf$search(
    query,
    sources = as.list(sources),
    query_type = query_type,
    limit = as.integer(limit)
  )

  .collection_to_tibble(result)
}


# Internal helper: convert a Python CompoundCollection to a tibble
.collection_to_tibble <- function(collection) {
  compounds <- collection$compounds
  if (length(compounds) == 0L) {
    return(tibble::tibble(
      cid = character(),
      name = character(),
      smiles = character(),
      inchikey = character(),
      formula = character(),
      molecular_weight = numeric()
    ))
  }

  rows <- lapply(compounds, function(c) {
    list(
      cid            = if (is.null(c$cid)) NA_character_ else as.character(c$cid),
      name           = if (is.null(c$name)) NA_character_ else as.character(c$name),
      smiles         = if (is.null(c$smiles)) NA_character_ else as.character(c$smiles),
      inchikey       = if (is.null(c$inchikey)) NA_character_ else as.character(c$inchikey),
      formula        = if (is.null(c$formula)) NA_character_ else as.character(c$formula),
      molecular_weight = if (is.null(c$properties) || is.null(c$properties$molecular_weight))
                           NA_real_
                         else as.numeric(c$properties$molecular_weight)
    )
  })

  tibble::as_tibble(do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE)))
}
