#' Compute molecular descriptors for one or more SMILES strings
#'
#' Uses RDKit (via the Python chemfuse package) to compute a comprehensive
#' set of 2D molecular descriptors including molecular weight, LogP, TPSA,
#' hydrogen bond donors/acceptors, and rotatable bonds.
#'
#' @param smiles Character vector. One or more SMILES strings.
#'
#' @return A tibble with one row per valid SMILES input and columns for each
#'   computed descriptor. Invalid SMILES result in rows with NA values.
#'
#' @examples
#' \dontrun{
#'   # Single compound
#'   desc <- cf_descriptors("CCO")
#'
#'   # Multiple compounds
#'   desc <- cf_descriptors(c("CCO", "CCC", "CC(=O)O"))
#' }
#'
#' @export
cf_descriptors <- function(smiles) {
  if (!is.character(smiles) || length(smiles) == 0L) {
    stop("'smiles' must be a non-empty character vector.")
  }

  desc_mod <- reticulate::import("chemfuse.compute.descriptors")

  rows <- lapply(smiles, function(smi) {
    tryCatch({
      result <- desc_mod$compute_descriptors(smi)
      py_dict <- reticulate::py_to_r(result)
      if (is.null(py_dict) || length(py_dict) == 0L) {
        list(smiles = smi, error = "invalid_smiles")
      } else {
        c(list(smiles = smi), py_dict)
      }
    }, error = function(e) {
      list(smiles = smi, error = conditionMessage(e))
    })
  })

  # Collect all descriptor names across all rows (to handle varying output)
  all_names <- unique(unlist(lapply(rows, names)))
  # Ensure consistent columns
  rows <- lapply(rows, function(r) {
    missing_cols <- setdiff(all_names, names(r))
    for (col in missing_cols) r[[col]] <- NA
    r[all_names]
  })

  result_df <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  tibble::as_tibble(result_df)
}


#' Predict ADMET properties for one or more SMILES strings
#'
#' Uses the ADMET-AI model (via the Python chemfuse package) to predict
#' absorption, distribution, metabolism, excretion, and toxicity properties.
#' Falls back to rule-based estimates if the ML model is unavailable.
#'
#' @param smiles Character vector. One or more SMILES strings.
#'
#' @return A tibble with one row per SMILES and columns for predicted
#'   ADMET properties including solubility, permeability, clearance,
#'   half-life, and toxicity endpoints.
#'
#' @examples
#' \dontrun{
#'   # Single compound
#'   admet <- cf_admet("CCO")
#'
#'   # Multiple compounds
#'   admet <- cf_admet(c("CCO", "CC(=O)Oc1ccccc1C(=O)O"))
#' }
#'
#' @export
cf_admet <- function(smiles) {
  if (!is.character(smiles) || length(smiles) == 0L) {
    stop("'smiles' must be a non-empty character vector.")
  }

  admet_mod <- reticulate::import("chemfuse.compute.admet")

  tryCatch({
    preds <- admet_mod$predict_admet(as.list(smiles))
    result <- reticulate::py_to_r(preds)
    if (is.data.frame(result)) {
      tibble::as_tibble(result)
    } else if (is.list(result)) {
      tibble::as_tibble(as.data.frame(result, stringsAsFactors = FALSE))
    } else {
      tibble::tibble(smiles = smiles, error = "unexpected_output")
    }
  }, error = function(e) {
    message("ADMET prediction failed: ", conditionMessage(e))
    tibble::tibble(smiles = smiles, error = conditionMessage(e))
  })
}


#' Compute drug-likeness filter results for a compound
#'
#' Evaluates a SMILES string against standard drug-likeness filters:
#' Lipinski Ro5, Veber, Ghose, Egan, and PAINS.
#'
#' @param smiles Character string. A single SMILES string.
#'
#' @return A named list with one element per filter, each containing
#'   a logical \code{passes} field and a character \code{reason} field.
#'
#' @examples
#' \dontrun{
#'   dl <- cf_druglikeness("CC(=O)Oc1ccccc1C(=O)O")
#'   dl$lipinski$passes  # TRUE
#' }
#'
#' @export
cf_druglikeness <- function(smiles) {
  if (!is.character(smiles) || length(smiles) != 1L || nchar(smiles) == 0L) {
    stop("'smiles' must be a single non-empty SMILES string.")
  }

  dl_mod <- reticulate::import("chemfuse.compute.druglikeness")

  tryCatch({
    result <- dl_mod$compute_all_filters(smiles)
    reticulate::py_to_r(result)
  }, error = function(e) {
    stop("Drug-likeness computation failed: ", conditionMessage(e))
  })
}
