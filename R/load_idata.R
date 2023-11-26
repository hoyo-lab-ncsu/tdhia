
#' load_idata
#'
#' \code{load_idata} returns the sum of all the values present in its arguments.
#'
#' This is a generic function: methods can be defined for it directly
#' or via the \code{\link{Summary}} group generic. For this to work properly,
#' the arguments \code{...} should be unnamed, and dispatch is on the
#' first argument.
#' @param idat_dir_paths the full filesystem path to directory containing the idat
#' files to be processed, formatted as a string or a vector of strings.
#' @param mft manifest file for particular imprintome array used in data collection.
#' A dataframe with metadata mapping probe_ids to CpG sites. See ?manifest_v1A2 dataset.
#' @param multicore boolean flag for whether multicore processing should be used
#'  when importing IDAT files.
#' @importFrom sesame openSesame
#' @importFrom stringr str_replace regex
#'
#'
#' @returns a matrix samples to load
#' @export
#'
#'
load_idata <- function(idat_dir_paths, mft = NULL, multicore = TRUE) {


  # Load manifest file
  mft = tdhia::manifest_v1A2
  print(mft)

  # Get list of IDAT files in target path
  all_idat_names <- dir(path = idat_dir_paths, pattern = "*.idat", ignore.case = TRUE)
  # Extract the basenames of IDAT files
  #   Removes _Grn.idat and _Red.idat from file names (case insensitive)
  all_idat_basenames <- stringr::str_replace(
    all_idat_names, pattern = stringr::regex("_(Grn|Red).idat$", ignore_case = TRUE),
    replacement = "")


  # Debugging: verify that two IDAT files are associated with each basename:
  #    One file for red and one for green channel fluorescence
  basename_tbl <- as.data.frame(table(all_idat_basenames))
  if (!all(basename_tbl$Freq==2)){
    warning(paste("IDAT_PAIR: Not all IDAT file basenames have 2 IDAT files associated with them.",
                  "Files are missing! See basename_tbl variable"))
  }

  # Get unique list of basenames (removes repeats)
  unq_idat_basenames <- unique(all_idat_basenames)
  #
  # # Load study meta_data
  # # study metadata file needs to have these columns
  # #   source: name of study that data originated from
  # #   patient_study_id: an integer, study-specific simple ID # for each patient
  # #   Patient_ID: the original medical record number/code for patient
  # #
  # # Note: the same patient participating in different studies will have different
  # #   Patient_IDs!
  # #
  # # Even though Patient_ID and patient_study_id have similar roles, its better to
  # #   Use a simple integer index for downstream processing.
  # # study_meta <- read_excel(study_meta_path)
  # # study_meta$patient_global_id <- paste(study_meta$source, study_meta$patient_study_id, sep = "_")
  #
  # # # List any patient IDs not found in IDAT files
  # # target_IDAT_basenames <- intersect(study_meta$Patient_ID, unq_idat_basenames)
  # # missing_idats <- setdiff(study_meta$Patient_ID, unq_idat_basenames)
  # # if (!length(missing_idats)==0) {
  # #   warning("There are missing IDAT files for some patients in study metadata, see 'missing_idats' variable")
  # # }
  #
  #
  #
  # # Processing IDAT files into beta matrices
  # ################################################################################
  # # Convert florescence data from array into beta values
  #
  # # Set Multicore parameter depending on operating system
  # if (.Platform$OS.type == "windows" && multicore) {
  #   multicore_arg <- BiocParallel::SnowParam(max(c(parallel::detectCores()-1,1)))
  # } else if (multicore) {
  #   multicore_arg <- BiocParallel::MulticoreParam(max(c(parallel::detectCores()-1,1)))
  # } else {
  #   multicore_arg = NULL
  # }
  #
  # # Load each of the IDAT file pairs, process with standard SeSame pipeline
  # probe_beta_matrix <- openSesame(paste0(idat_dir_paths, '/', target_IDAT_basenames),
  #                                 platform="TruDx_imprintome", manifest = mft, mask=F,
  #                                 BPPARAM = multicore_arg)
  #
  # return(probe_beta_matrix)

}

