
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
#' @param idat_dir_paths input idat files that can be in two forms: (1) the full
#' filesystem path to directory containing the idat files to be processed,
#' formatted as a string or a vector of strings. (2) a vector of sigset objects
#' (idat files once they are imported into memory, this option is used for
#' simualted idat files since we don't save them to disk).
#' @param platform string for the platform that the array belongs to, as specified
#' within the SeSame package with openSesame(). default: "TruDx_imprintome"
#' @param mft manifest file for particular imprintome array used in data collection.
#' A dataframe with metadata mapping probe_ids to CpG sites and genome locations.
#' See ?manifest_v1A2 for more info.
#' @param multicore boolean flag or integer, when true the max number of cores minus 1 is
#' used for the current system, FALSE sets to single core operation, any integer
#' greater than zero sets to the specified number of cores.
#' @param sesame_prep string of number/ letters to control sesame normalization
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
#' @param quantile_norm a boolean flag, when set to TRUE, applies quantile
#' normalization between the columns in the probe_beta dataframe (default = FALSE).
#' @param mask a boolean flag specified in opensesame() in the Sesame package,
#' when TRUE excludes some probes due to issues of inter-dependence of measurements
#' (default = FALSE).
#' @param db_flag boolean when true exports function workspace to disk.
#' @param enforce_req_idats check that all idat files requested in idat_basenames
#'  are found on disk. Throws error if this is not the case.
#' @param enforce_idat_names boolean when true raises error if input IDAT files
#' do not follow proper capitalization pattern, if FALSE only issues warning.
#' @returns a named list with the following fields:
#'  - probe_beta_df: dataframe of beta values, probe_id x sample_id
#'  - probe_pval_df: dataframe of signal p-values, probe_id x sample_id
#'  - platform: string that describes platform for methylation array
#'  - manifest: dataframe of the manifest file used for SeSame processing
#' @export
load_idata_to_probes <-
  function(idat_dir_paths, platform = "TruDx_imprintome",
           mft = NULL, multicore = TRUE, sesame_prep = "0CDB",
           idat_basenames = NULL, quantile_norm = FALSE, mask = FALSE,
           db_flag = FALSE, enforce_req_idats = FALSE, enforce_idat_names = TRUE) {
    # Load manifest file if platform is true diagnostic imprintome array
    if (base::is.null(mft) && platform=="TruDx_imprintome") {mft = tdhia::manifest_v1A2}

    if(db_flag) save(list = ls(all.names = TRUE), file = "load_idata_to_probes.RData")
    # load(file = "load_idata_to_probes.RData")


    if (!methods::is(idat_dir_paths[[1]],"data.frame")) {
      # Get list of IDAT files in target path
      obs_idat_fullnames <- dir(path = idat_dir_paths, pattern = "/*.idat", full.names = TRUE, ignore.case = TRUE)

      # verify that IDAT files are properly named, raise error if not
      correct_idat_names(obs_idat_fullnames, rename = FALSE, failcheck_error = enforce_idat_names)


      # Find full path and basename of observed IDAT files
      #   Removes _Grn.idat and _Red.idat from file names (case insensitive)
      obs_idat_full_basenames <- stringr::str_replace(
        obs_idat_fullnames, pattern = stringr::regex("_(Grn|Red).idat$", ignore_case = TRUE),
        replacement = "")
      obs_idat_basenames <- basename(obs_idat_full_basenames)


      # Debugging: verify that two IDAT files are associated with each basename:
      #    One file for red and one for green channel fluorescence
      fullname_tbl <- as.data.frame(table(obs_idat_full_basenames))
      fullname_tbl$obs_idat_full_basenames <- as.character(fullname_tbl$obs_idat_full_basenames)

      basename_tbl <- as.data.frame(table(obs_idat_basenames))
      basename_tbl$obs_idat_basenames <- as.character(basename_tbl$obs_idat_basenames)

      # Verify and error if not all requested idats are included, but could be
      # split across dirs
      if (enforce_req_idats) {
        is_requested <- basename(basename_tbl$obs_idat_basenames) %in%
          idat_basenames
        if (sum(is_requested) != length(idat_basenames)) {
          stop("Not all requested idat files were found in directories.")
        }
      }

      if (!all(basename_tbl$Freq==2)) {
        warning(paste("IDAT_PAIR: Not all IDAT file basenames have 2 IDAT files associated with them.",
                      "Files are missing! See basename_tbl variable"))
      }

      # Get unique list of basenames (removes repeats) found in input folder
      unq_obs_idat_basenames <- unique(obs_idat_full_basenames)

      # Check that all specified idat_basenames specified exist in unq_obs_idat_basenames
      if (!all(idat_basenames %in% basename(unq_obs_idat_basenames)) & enforce_req_idats) {
        stop("load_idata_to_probes: some idat_basenames not found.")
      }

      # Discard IDAT files that are not included in idat_basenames list
      if (!is.null(idat_basenames)) {
        if (is.data.frame(idat_basenames)) {
          sub_idats <- idat_basenames$idat_basename
        } else if (is.vector(idat_basenames)) {
          sub_idats <- idat_basenames
        } else {stop ("idat_basenames need to be a vector or dataframe")}

        is_requested <- is.element(basename(unq_obs_idat_basenames), sub_idats)

        unq_obs_idat_basenames <- unq_obs_idat_basenames[is_requested]
      }

    } else {
      unq_obs_idat_basenames <- idat_dir_paths
    }


    # SeSame Processing
    #___________________________________________________________________________

    # Caching SeSame Data Files (required for running SeSame firs time...)
    cat("Checking that Sesame idatSignature file is locally cached...")
    sesameData::sesameDataCache()
    cat("Done.\n")

    # Set Multicore parameter depending on operating system
    if (multicore > 1) {
      default_ncores <- multicore
    } else {
      default_ncores <- max(c(parallel::detectCores()-1,1))
    }
    if (.Platform$OS.type == "windows" && multicore) {
      multicore_arg <- BiocParallel::SnowParam(default_ncores)
    } else if (multicore) {
      multicore_arg <- BiocParallel::MulticoreParam(default_ncores,progressbar = TRUE)
    } else {
      multicore_arg = NULL
    }


    # Load each of the IDAT file pairs, process with standard SeSame pipeline
    cat(sprintf("Processing %.0f IDAT Files...", length(unq_obs_idat_basenames)))
    if (!is.null(multicore_arg)) {
      cat(sprintf(" (Using %i/%i cores.)...", default_ncores, parallel::detectCores()))
    }

    idat.out <- process_IDATS(unq_obs_idat_basenames, platform, mft, multicore_arg)
    probe_beta_matrix = idat.out$betas
    probe_pval_matrix = idat.out$pvals


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


    probe_beta <- list(probe_beta_df = probe_beta_df,
                       probe_pval_df = probe_pval_df,
                       platform = platform,
                       manifest = mft)
    return(probe_beta)

  }



#' process_idats
#'
#' Processes idat files
#'
#' beta = M/(M+U)
#'
#' This function processes the signal for all probes in the array, with no
#' filtering performed.
#'
#' @param unq_obs_idat_basenames list of idat basenames, or list of imported
#' sigsets, which is SeSame's class for imported IDAT files
#' @param platform string for the platform that the array belongs to, as specified
#' within the SeSame package with openSesame(). default: "TruDx_imprintome"
#' @param mft manifest file for particular imprintome array used in data collection.
#' A dataframe with metadata mapping probe_ids to CpG sites and genome locations.
#' See ?manifest_v1A2 for more info.
#' @param multicore_arg boolean flag, when true the max number of cores minus 1 is
#' used for the current system.
#'
#' @returns a named list with the following fields:
#'  - betas: matrix of methylation beta values, probes x patients
#'  - pvals: matrix of signal p-values, probes x patients
#'
process_IDATS <- function(unq_obs_idat_basenames, platform, mft, multicore_arg) {


  if (is.character(unq_obs_idat_basenames)) {
    ss = BiocParallel::bplapply(
      unq_obs_idat_basenames, function(pfx) {
        sesame::readIDATpair(pfx,  platform = "", manifest = mft)}, BPPARAM = multicore_arg)
  } else {
    ss = unq_obs_idat_basenames
  }

  betas = do.call(cbind,BiocParallel::bplapply(ss, function(ss) {
    sesame::getBetas(
      sdf =  sesame::noob(
        sdf = sesame::dyeBiasNL(
          sdf = sesame::inferInfiniumIChannel(ss))))},
    BPPARAM = multicore_arg))

  pvals <- do.call(cbind,BiocParallel::bplapply(ss, function(ss) {
    sesame::pOOBAH(return.pval = TRUE,
                   sdf =  sesame::dyeBiasNL(
                     sdf =  sesame::inferInfiniumIChannel(ss)))},
    BPPARAM = multicore_arg))

  return(list(betas = betas, pvals = pvals))
}
