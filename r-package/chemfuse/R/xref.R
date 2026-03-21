#' Cross-database identifier mapping
#'
#' Maps a compound identifier (PubChem CID, ChEMBL ID, SMILES, or InChIKey)
#' to identifiers in other databases using the UniChem service.
#'
#' @param cid Integer or character. PubChem CID. Optional.
#' @param chembl_id Character string. ChEMBL identifier (e.g.,
#'   \code{"CHEMBL25"}). Optional.
#' @param smiles Character string. SMILES string. Optional.
#' @param inchikey Character string. Standard InChIKey (27 chars). Optional.
#'
#' At least one of \code{cid}, \code{chembl_id}, \code{smiles}, or
#' \code{inchikey} must be provided.
#'
#' @return A tibble with one row per database mapping. Columns include
#'   \code{database} (database name) and \code{identifier} (mapped ID).
#'   Common databases include \code{pubchem}, \code{chembl},
#'   \code{chebi}, \code{drugbank}, \code{kegg}, \code{hmdb}.
#'
#' @examples
#' \dontrun{
#'   # Map by PubChem CID (aspirin = 2244)
#'   xref <- cf_xref(cid = 2244)
#'
#'   # Map by ChEMBL ID
#'   xref <- cf_xref(chembl_id = "CHEMBL25")
#'
#'   # Map by SMILES
#'   xref <- cf_xref(smiles = "CC(=O)Oc1ccccc1C(=O)O")
#' }
#'
#' @export
cf_xref <- function(
    cid = NULL,
    chembl_id = NULL,
    smiles = NULL,
    inchikey = NULL
) {
  # Validate at least one argument provided
  if (is.null(cid) && is.null(chembl_id) && is.null(smiles) && is.null(inchikey)) {
    stop("At least one of 'cid', 'chembl_id', 'smiles', or 'inchikey' must be provided.")
  }

  # Build keyword arguments for Python call
  kwargs <- list()
  if (!is.null(cid))       kwargs$cid       <- as.integer(cid)
  if (!is.null(chembl_id)) kwargs$chembl_id  <- as.character(chembl_id)
  if (!is.null(smiles))    kwargs$smiles     <- as.character(smiles)
  if (!is.null(inchikey))  kwargs$inchikey   <- as.character(inchikey)

  result <- do.call(cf$map_identifiers, kwargs)
  mapping <- reticulate::py_to_r(result)

  if (length(mapping) == 0L) {
    return(tibble::tibble(database = character(), identifier = character()))
  }

  tibble::tibble(
    database   = names(mapping),
    identifier = unlist(mapping, use.names = FALSE)
  )
}
