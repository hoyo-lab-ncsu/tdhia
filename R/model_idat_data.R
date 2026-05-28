

#' model_idat_data
#'
#' Produces a rough model to simulate IDAT data from a reference IDAT file.
#' Can be used to test pipelines or vignettes for the package. Simulated idat
#' files can be produced with the returned data object and calling:
#'
#' simulate_idat_data().
#'
#' Data is stored as a simple ecfd, later versions should switch to a smoothed
#' pdf of the fluorescent values.
#'
#' @param ref_idat_path full path to idat file, including prefix of idat.
#' @param model_save_path if specified, saves model to rda file.
#' @param data_colnames  column names for fluorescent data
#' @return sigset_model, a named list with the follow fields:
#'    sigset_template: a data frame that is copied from the reference idat data,
#'       with all the numerical data cleared.
#'    sigset_mask: a boolean mask specifying what rows of data columns were NA.
#'    sigset_ecdfs: a set of empirical distribution functions for each data
#'       column.
#' @export
model_idat_data <- function(ref_idat_path, model_save_path = NA, data_colnames =
                              c("MG", "MR", "UG", "UR")) {

# Load IDAT data
ref_sigset <- sesame::readIDATpair(ref_idat_path, platform = "TDHIA",
  manifest = tdhia::manifest_v1A2, controls = NULL, verbose = FALSE)


# We will store modeling data about sigset in a named list.
sigset_model <- list()
# Keep the Probe_ID, col, mask, and attributes directly.
sigset_model$sigset_template <- ref_sigset
sigset_model$sigset_template[,data_colnames] <- NA

# Save mask of what values are NA for data
sigset_model$sigset_mask <- as.data.frame(
  is.na(dplyr::select(ref_sigset, dplyr::any_of(data_colnames))))


# Create an empirical distribution function for each column
sigset_ecdfs = list()
for (n in seq_along(data_colnames)) {
  sigset_model$sigset_values[[data_colnames[n]]] =
    sort(na.omit(ref_sigset[,data_colnames[n]]))
}

 if (!is.na(model_save_path)) {
   save(sigset_model, file = paste0(model_save_path, "/sigset_model.rda"))
 }

  return(sigset_model)
}
