
#' load_idata_to_probes
#'
#' Loads a series flourescence measurements from a custom infinium methylation array
#' and caclulates the methylation beta values and associated detection p-values.
#'
#' beta = M/(M+U)
#'
#' This function processes the signal for all probes in the array, with no
#' filtering performed.
#'
#' @param idat_dir_paths the full filesystem path to directory containing the idat
#' files to be processed, formatted as a string or a vector of strings.
#'
#' @param platform string for the platform that the array belongs to, as specified
#' within the SeSame package with openSesame(). default: "TruDx_imprintome"
#'
#' @param mft manifest file for particular imprintome array used in data collection.
#' A dataframe with metadata mapping probe_ids to CpG sites and genome locations.
#' See ?manifest_v1A2 for more info.
#'
#' @param multicore boolean flag, when true the max number of cores minus 1 is
#' used for the current system.
#'
#' @param sesame_pre string of number/ letters to control sesame normalization
#' steps. Default is 0CDB.
#' @param idat_basenames a character vector or dataframe.
#' If idat_basenames is a character vector, then the elements are the
#' idat_basenames to be processed. An idat basename is the filename without the
#' red and green channel suffix included. So the idat basename
#' "207344530007_R08C02" would be used for the pair of files:
#'   1. 207344530007_R08C02_Red.idat
#'   2. 207344530007_R08C02_Grn.idat
#' If idat_basenames is a dataframe, then it contains two columns: idat_basename
#' and id. idat_basename is a character vector of idat basenames. id is a numeric
#'  or string vector for what each idat basename should be renamed to.
#'
#' @param quantile_norm a boolean flag, when set to TRUE, applies quantile
#' normalization between the columns in the probe_beta dataframe (default = FALSE).
#'
#' @param mask a boolean flag specified in opensesame() in the Sesame package,
#' when TRUE excludes some probes due to issues of inter-dependence of measurements
#' (default = FALSE).
#'
#' @returns a named list with the following fields:
#'  - probe_beta_df: dataframe of beta values, probe_id x sample_id
#'  - probe_pval_df: dataframe of signal p-values, probe_id x sample_id
#'  - platform: string that describes platform for methylation array
#'  - manifest: dataframe of the manifest file used for SeSame processing
#'
load_idata_to_probes <-
  function(idat_dir_paths, platform = "TruDx_imprintome",
           mft = NULL, multicore = TRUE, sesame_prep = "0CDB",
           idat_basenames = NULL, quantile_norm = FALSE, mask = FALSE) {
  # Load manifest file if platform is true diagnostic imprintome array
  if (is.null(mft) && platform=="TruDx_imprintome") {mft = tdhia::manifest_v1A2}
  # save(list = ls(all.names = TRUE), file = "load_idata_debug.RData")
  # load(file = "load_idata_debug.RData")

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

  # Get unique list of basenames (removes repeats) found in input folder
  obs_idat_basenames <- unique(all_idat_basenames)

  # Discard IDAT files that are not included in idat_basenames list
  if (!is.null(idat_basenames)) {
    if (is.data.frame(idat_basenames)) {
      sub_idats <- idat_basenames$idat_basename
    } else if (is.vector(idat_basenames)) {
      sub_idats <- idat_basenames
    } else {stop ("idat_basenames need to be a vector or dataframe")}
    obs_idat_basenames <-
      obs_idat_basenames[is.element(obs_idat_basenames, sub_idats)]
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
  cat(sprintf("Processing %.0f IDAT Files...", length(obs_idat_basenames)))
  if (!is.null(multicore_arg)) {
    cat(sprintf(" (Using %i/%i cores.)...", default_ncores, parallel::detectCores()))
  }
  probe_beta_matrix <-
    sesame::openSesame(paste0(idat_dir_paths, '/', obs_idat_basenames),
                       platform=platform, manifest = mft, prep = sesame_prep,
                       BPPARAM = multicore_arg, fun = sesame::getBetas)


  probe_pval_matrix <-
    sesame::openSesame(paste0(idat_dir_paths, '/', obs_idat_basenames),
                       platform=platform, manifest = mft, prep = sesame_prep,
                       BPPARAM = multicore_arg, fun = sesame::pOOBAH,
                       return.pval  =TRUE)

  # Debugging: saving output from sesame pipeline
  # Todo: pipeline is currently run twice to get beta and p-values, recode to run once
  # save(list = ls(all.names = TRUE), file = "load_idata_debug.RData")
  # load(file = "load_idata_debug.RData")


  # Convert matrix to dataframe for easier manipulation
  probe_beta_df <- as.data.frame(probe_beta_matrix)
  probe_pval_df <- as.data.frame(probe_pval_matrix)
  cat("Done.\n")

  # If a column mapping is specified, rename column names in probe_beta_df
  #   Discard unmapped columns if user specifies it
  if (is.data.frame(idat_basenames)) {

      probe_beta_df <- probe_beta_df[,is.element(colnames(probe_beta_df),
                                  idat_basenames$idat_basename)]
      probe_pval_df <- probe_pval_df[,is.element(colnames(probe_pval_df),
                                                 idat_basenames$idat_basename)]

    # Rename by matching new id to each idat basename using a left_join()
    remapped_colnames <-
      dplyr::left_join(df_remap <- data.frame(idat_basename = colnames(probe_beta_df)),
                       dplyr::select(idat_basenames, c("idat_basename", "id")),
                       by = "idat_basename",
                     na_matches = "never", unmatched = "error",
                     relationship = "one-to-one")
    # Rename columns of probe_beta_df and probe_pval_df
    colnames(probe_beta_df) <- remapped_colnames$id
    colnames(probe_pval_df) <- remapped_colnames$id
  }


  # Apply quantile normalization between all columns if specified by user.
  if (quantile_norm) {
    norm_probe_df <-
      as.data.frame(preprocessCore::normalize.quantiles(as.matrix(probe_beta_df)),
                    row.names = rownames(probe_beta_df))
    colnames(norm_probe_df) <- colnames(probe_beta_df)
    probe_beta_df <- norm_probe_df
  }

  # # Add probe_id as a first column and remove row names
  # probe_beta_df<-cbind(data.frame(Probe_ID = rownames(probe_beta_df)),probe_beta_df)
  # rownames(probe_beta_df) <- NULL
  # probe_pval_df<-cbind(data.frame(Probe_ID = rownames(probe_pval_df)),probe_pval_df)
  # rownames(probe_pval_df) <- NULL



  probe_beta <- list(probe_beta_df = probe_beta_df,
                      probe_pval_df = probe_pval_df,
                      platform = platform,
                      manifest = mft)
  return(probe_beta)

}


