
#' load_idata
#'
#' Returns the sum of all the values present in its arguments.
#'
#' @param idat_dir_paths the full filesystem path to directory containing the idat
#' files to be processed, formatted as a string or a vector of strings.
#' @param mft manifest file for particular imprintome array used in data collection.
#' A dataframe with metadata mapping probe_ids to CpG sites. See ?manifest_v1A2 dataset.
#' @param multicore boolean flag for whether multicore processing should be used
#'  when importing IDAT files.
#' @param column_mapping a dataframe containing two columns used to convert the
#' basename of the IDAT files (used to label columns of output dataframe), into
#' a different set of identifiers/mappings.
#'  - idat_basename: column of base filenames of IDAT (no repeated values).
#'  - id: column of strings of new mapping (no repeated values).
#' @param discard_unmapped a boolean flag when set to true discards any idat base
#'  filenames that are not included in the column_mapping argument (note:
#'  column_mapping must be specified).
#' @param quantile_norm a boolean flag, when set to TRUE, applies quantile
#' normalization between the columns in the returned dataframe.
#'
#' @returns a dataframe containing beta values where rows are probe_ids and columns
#' are the basenames of the idat files or some other mapping specifed by column_mappings.
#'
#' @export
#'
load_idata <- function(idat_dir_paths, mft = NULL, multicore = TRUE,
                       column_mapping = NULL, discard_unmapped = TRUE,
                       quantile_norm = TRUE) {


  # Load manifest file
  if (is.null(mft)) {mft = tdhia::manifest_v1A2}

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
  idat_basenames <- unique(all_idat_basenames)

  # Discard IDAT files that are not included in column_mapping file
  if (!is.null(column_mapping) &&  discard_unmapped) {
    idat_basenames <- idat_basenames[is.element(idat_basenames,
                                                column_mapping$idat_basename)]
  }


  # Caching SeSame Data Files (required for running SeSame firs time...)
  cat("Checking that Sesame idatSignature file is locally cached...")
  sesameData::sesameDataCache()
  cat("Done.\n")

  # Set Multicore parameter depending on operating system
  default_ncores <- max(c(parallel::detectCores()-1,1))
  if (.Platform$OS.type == "windows" && multicore) {
    multicore_arg <- BiocParallel::SnowParam(default_ncores)
  } else if (multicore) {
    multicore_arg <- BiocParallel::MulticoreParam(default_ncores)
  } else {
    multicore_arg = NULL
  }

  # Load each of the IDAT file pairs, process with standard SeSame pipeline
  cat("Processing IDAT Files...")
  if (!is.null(multicore_arg)) {
    cat(sprintf(" (Using %i/%i cores.)...", default_ncores, parallel::detectCores()))
  }
  probe_beta_matrix <- sesame::openSesame(paste0(idat_dir_paths, '/', idat_basenames),
                                  platform="TruDx_imprintome", manifest = mft, mask=F,
                                  BPPARAM = multicore_arg)

  # Convert matrix to dataframe for easier manipulation
  probe_beta_df <- as.data.frame(probe_beta_matrix)
  cat("Done.\n")

  # If a column mapping is specified, rename column names in probe_beta_df
  #   Discard unmapped columns if user specifies it
  if (!is.null(column_mapping)) {

    # If discard_unmapped is TRUE, any columns in [probe_beta_df] not included in
    # [column mapping] are discarded. This is useful if the IDAT folder has data
    # from multiple studies
    if (discard_unmapped) {
      probe_beta_df <- probe_beta_df[,is.element(colnames(probe_beta_df),
                                  column_mapping$idat_basename)]
    }

    # Match new id to each idat_basename using a left join
    remapped_colnames <-
      dplyr::left_join(df_remap <- data.frame(idat_basename = colnames(probe_beta_df)),
                       dplyr::select(column_mapping, c("idat_basename", "id")),
                       by = "idat_basename",
                     na_matches = "never", unmatched = "error", relationship = "one-to-one")
    # Rename columns of probe_beta_df
    colnames(probe_beta_df) <- remapped_colnames$id
  }

  # Apply quantile normalization between all columns if specified by user.
  if (quantile_norm) {
    norm_probe_df <-
      as.data.frame(preprocessCore::normalize.quantiles(as.matrix(probe_beta_df)),
                    row.names = rownames(probe_beta_df))
    colnames(norm_probe_df) <- colnames(probe_beta_df)
    prob_beta_df <- norm_probe_df
  }

  return(probe_beta_df)

}


