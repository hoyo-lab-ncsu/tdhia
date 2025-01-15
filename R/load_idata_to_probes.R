
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
#' @param merge_replicates string with the following possible values that
#' determines how replicate probes are handled.
#'  - NULL: the probe replicates are not merged (default value).
#'  - pre_beta: replicates are merged by averaging fluorescent signal from
#'   each channel individually, beta values are then calculated from these
#'   averaged values (pre_beta = merging done before beta calculation).
#'  - post_beta: beta value is calculated before replicates are merged and beta
#'   values are averaged between replicates (post_beta = merging done before
#'   beta calculation).
#' @returns a named list with the following fields:
#'  - probe_beta_df: dataframe of beta values, probe_id x sample_id
#'  - probe_pval_df: dataframe of signal p-values, probe_id x sample_id
#'  - platform: string that describes platform for methylation array
#'  - manifest: dataframe of the manifest file used for SeSame processing
#' @export
load_idata_to_probes <-
  function(idat_dir_paths, platform = "TruDx_imprintome",
           mft = NULL, multicore = 1,
           idat_basenames = NULL, quantile_norm = FALSE, mask = FALSE,
           db_flag = FALSE, enforce_req_idats = FALSE, enforce_idat_names = TRUE,
           merge_replicates = NULL) {
    # Load manifest file if platform is true diagnostic imprintome array
    if (base::is.null(mft) && platform=="TruDx_imprintome") {mft = tdhia::manifest_v1A2_design_scores}

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

    # Set multicore parameter depending on operating system
    if (is.logical(multicore) && multicore) {
      multicore <- max(c(parallel::detectCores()-1,1))
    }

    if (!multicore || multicore==1) {
      core_params <- BiocParallel::SerialParam()
    } else if (.Platform$OS.type == "windows") {
      core_params <- BiocParallel::SnowParam(multicore,)
    }  else {
      core_params <- BiocParallel::MulticoreParam(multicore, progressbar = TRUE)
    }


    # Load each of the IDAT file pairs, process with standard SeSame pipeline
    cat(sprintf("Processing %.0f IDAT Files...", length(unq_obs_idat_basenames)))
    if (!is.null(core_params)) {
      cat(sprintf(" (Using %i/%i cores.)...", core_params$workers, parallel::detectCores()))
    }

    idat.out <- process_IDATS(unq_obs_idat_basenames, platform = platform,
                              mft = mft, core_params = core_params,
                              db_flag = db_flag, merge_replicates = merge_replicates)
    probe_beta_matrix = idat.out$betas
    probe_pval_matrix = idat.out$p_vals


    # if (!is.character(unq_obs_idat_basenames)) {
    #   sample_names <- paste0("Sample_", length(unq_obs_idat_basenames))
    # } else {
    #   sample_names <- basename(unq_obs_idat_basenames)
    # }
    # colnames(probe_beta_matrix) <- sample_names
    # colnames(probe_pval_matrix) <- sample_names

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
#' @description
#' Processes idat files
#'
#' beta = M/(M+U)
#'
#' This function processes the signal for all probes in the array, with no
#' filtering performed.
#' Sesame prep code: 0CDB
#' @param unq_obs_idat_basenames list of idat basenames, or list of imported
#' sigsets, which is SeSame's class for imported IDAT files
#' @param platform string for the platform that the array belongs to, as specified
#' within the SeSame package with openSesame(). default: "TruDx_imprintome"
#' @param mft manifest file for particular imprintome array used in data collection.
#' A dataframe with metadata mapping probe_ids to CpG sites and genome locations.
#' See ?manifest_v1A2 for more info.
#' @param merge_replicates string with the following possible values that
#' determines how replicate probes are handled.
#'  - NULL: the probe replicates are not merged (default value).
#'  - pre_beta: replicates are merged by averaging fluorescent signal from
#'   each channel individually, beta values are then calculated from these
#'   averaged values (pre_beta = merging done before beta calculation).
#'  - post_beta: beta value is calculated before replicates are merged and beta
#'   values are averaged between replicates (post_beta = merging done before
#'   beta calculation).
#' @param core_params boolean flag, when true the max number of cores minus 1 is
#' used for the current system.
#' @param db_flag todo
#' @returns a named list with the following fields:
#'  - betas: matrix of methylation beta values, probes x patients
#'  - p_vals: matrix of signal p-values, probes x patients
#'
process_IDATS <- function(unq_obs_idat_basenames, platform, mft, core_params,
                          merge_replicates = NULL, db_flag = FALSE) {

  if(db_flag) save(list = ls(all.names = TRUE), file = "process_IDATS.RData")
  # load(file = "process_IDATS.RData")


  if (is.character(unq_obs_idat_basenames)) {
    ss = BiocParallel::bplapply(
      unq_obs_idat_basenames, function(pfx) {
        sesame::readIDATpair(pfx,  platform = "", manifest = mft)},
      BPPARAM = BiocParallel::SerialParam())
    names(ss) <- basename(unq_obs_idat_basenames)
  } else {
    ss = unq_obs_idat_basenames
    unq_obs_idat_basenames = paste0("Sample_", 1:length(unq_obs_idat_basenames))
  }

  # Calculate sesame p-values for individual probes
  p_vals <- do.call(cbind, BiocParallel::bplapply(ss, function(ss) {
    sesame::pOOBAH(return.pval = TRUE,
                   sdf =  sesame::dyeBiasNL(
                     sdf =  sesame::inferInfiniumIChannel(ss)))},
    BPPARAM = core_params))

  # colnames(p_vals) <- basename(unq_obs_idat_basenames)

  # Perform dye bias corection and background subtraction
  ss_sig <- BiocParallel::bplapply(ss, function(ss) {
    ss_sig =  sesame::noob(
      sdf = sesame::dyeBiasNL(
        sdf = sesame::inferInfiniumIChannel(ss)))},
    BPPARAM = core_params)


  ## Add p_vals to sigset list
  ss_sig <- lapply(1:length(ss_sig), function(x) ss_sig[[x]] %>%
                     dplyr::mutate(p_val = p_vals[,x]))
  names(ss_sig) <- names(ss)


  if(db_flag) save(list = ls(all.names = TRUE), file = "process_IDATS.RData")
  # load(file = "process_IDATS.RData")

  if (!is.null(merge_replicates)) {
    # Get list of non-unique CPG IDs to be merged
    # Calculate merged beta (pre and post merge) for
    # merged_data <- BiocParallel::bplapply(
    #   ss_sig, function(x)
    #     sigset_merge_replicates(x, merge_replicates = merge_replicates),
    #   BPPARAM = BiocParallel::SerialParam())

    merged_data_list = lapply(ss_sig, function(x)
      sigset_merge_replicates(x, merge_replicates = merge_replicates, db_flag = db_flag))

    # Calculate beta values from merged sig sets
    betas = do.call(cbind, lapply(merged_data_list, FUN = function(x) x$betas))
    p_vals = do.call(cbind, lapply(merged_data_list, FUN = function(x) x$p_vals))

    #
    # colnames(p_vals) <- basename(unq_obs_idat_basenames)

  } else {

    # Calculate beta values for all probes (no merging or beta and p_val with probes)
    betas = do.call(cbind,BiocParallel::bplapply(ss_sig, function(ss) {
      sesame::getBetas(ss)},
      BPPARAM = BiocParallel::SerialParam()))
    colnames(betas) <- basename(unq_obs_idat_basenames)
  }

  out <- list(betas = betas, p_vals = p_vals)

  return(out)
}



#' sigset_merge_replicates
#' @description
#' merges the signal and p-values of probe replicates. signal is averaged by
#' fluorescent channel individually, p-values are merged based on a bates distribution.
#'
#' @param sigset dataframe that is a sigset from the sesame package.
#' @param merge_replicates string with the following possible values that
#' determines how replicate probes are handled.
#'  - NULL: the probe replicates are not merged (default value).
#'  - pre_beta: replicates are merged by averaging fluorescent signal from
#'   each channel individually, beta values are then calculated from these
#'   averaged values (pre_beta = merging done before beta calculation).
#'  - post_beta: beta value is calculated before replicates are merged and beta
#'   values are averaged between replicates (post_beta = merging done before
#'   beta calculation).
#' @param db_flag boolean, when TRUE the workspace is saved to disk for debugging.
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#' @export
sigset_merge_replicates <- function(sigset, merge_replicates = "post_merge", db_flag = TRUE) {

  if(db_flag) save(list = ls(all.names = TRUE), file = "sigset_merge_replicates.RData")
  # load(file = "sigset_merge_replicates.RData")

  if (!merge_replicates %in% c(NULL, "pre_beta", "post_beta")) {
    stop("argument merge_replicates needs to be NULL, 'pre_beta','or 'post_beta'.")
  }

  # Initialize variables to
  cpg_ids <- Probe_ID <- MR <- MG <- UG <- UR <- mask <- p_val <- key <- NA

  # Add cpg_ids and key column (for sorting rows below to preserve order)
  probe_ids <- sigset$Probe_ID
  sigset$cpg_ids <- stringr::str_replace(string = probe_ids, pattern = "_.{4}$", replacement = "")
  sigset$key = 1:nrow(sigset)
  # cpg_table <- table(sigset$cpg_ids)
  # nunq_cpg <- names(cpg_table)[cpg_table > 1]

  # Calculate beta signal before any probe replicate averaging
  sigset$beta_orig <- sesame::getBetas(sigset)

  # Probe renaming function for merged replicates
  fProbe_ID <- function(x) {if (length(x)==1) {y=x[1]
  } else {y=stringr::str_replace(x[1], "_.{4}$", "_MERG")}; return(y)}

  # Calculate mean fluorescent signal across replicates, record original values also
  sigset_merge <- sigset %>% dplyr::group_by(cpg_ids) %>% dplyr::summarize(
    # Rename this column later, prevent wrong sub in
    Probe_ID_ = fProbe_ID(.data$Probe_ID),
    MG = mean(MG),  MR = mean(MR), UG = mean(UG), UR = mean(UR),
    col = col[1],  mask = mask[1], mean_orig_beta = mean(.data$beta_orig),
    mean_p_val = mean(p_val),
    # Probe_ID_ renamed above to ensure this statement has correct scope
    orig_probe_ids = paste(.data$Probe_ID,collapse=" "),
    orig_p_vals = paste(p_val,collapse=" "),
    orig_betas = paste(beta_orig,collapse=" "),
    n = length(MG), key = min(key),
    n = length(p_val)) %>%
    dplyr::arrange(key) %>% dplyr::rename(Probe_ID = .data$Probe_ID_)

  # Calculate beta values based on merged fluorescent values
  sigset_merge$beta_merged <- sesame::getBetas(sigset_merge)

  # Calculate p-values of the merged signal between probes
  sigset_merge$merged_p_val = pbates(mean_p_vals = sigset_merge$mean_p_val, samples = sigset_merge$n)

  # Fill in default values for final probe_id, beta, p_val, and n fields
  if (merge_replicates == "pre_merge") {
    sigset_merge$final_beta <- sigset_merge$mean_orig_beta
  } else {sigset_merge$final_beta <- sigset_merge$beta_merged}
  sigset_merge$final_Probe_ID <- sigset_merge$Probe_ID
  sigset_merge$final_p_val <-  sigset_merge$merged_p_val
  sigset_merge$final_n <- sigset_merge$n

  # Scan for merged probes
  rep_ind <- which(sigset_merge$n > 1)

  # For the merged probes select best measurement (merged or individual)
  for (k in seq_along(rep_ind)) {
    # For each row, assemble vectors of merged and individual betas,  p_vals, and probe_ids
    rind = rep_ind[k]
    # Assemble vector of merged and individual betas
    if (merge_replicates == "pre_merge") {merged_beta <- sigset_merge$mean_orig_beta[rind]
    } else {merged_beta <- sigset_merge$beta_merged[rind]}
    betas <- c(merged_beta, sapply(stringr::str_split(
      sigset_merge$orig_betas[rind], pattern = " "), as.double))

    # Assemble vector of merged and individual p_values
    p_vals <- c(sigset_merge$merged_p_val[rind], sapply(stringr::str_split(
      sigset_merge$orig_p_vals[rind], pattern = " "), as.double))

    # Assemble vector of merged and individual probe_ids
    probe_ids <-  c(sigset_merge$Probe_ID[rind], stringr::str_split(
      sigset_merge$orig_probe_ids[rind], pattern = " ")[[1]])

    # Determine if merged probe or individual has best p_value
    ind <- which.min(p_vals)
    if (length(ind)==0) ind = length(p_vals)
    # Select betas, p_vals, and probe_id using this index


    # Update final probe_id
    sigset_merge$final_Probe_ID[rind] <- probe_ids[ind]
    # Update final betas
    sigset_merge$final_beta[rind] <- betas[ind]
    # Update final p value
    sigset_merge$final_p_val[rind] <- p_vals[ind]
    # If one probe is used, change number of probes to 1
    if (ind!=1) sigset_merge$final_n[rind] <- 1
  }

  # Add probe_ids as rownames to beta and p_val values
  names(sigset_merge$final_beta) <- sigset_merge$final_Probe_ID
  names(sigset_merge$final_p_val) <- sigset_merge$final_Probe_ID
  # Export beta and p_val
  out <- list(betas = sigset_merge$final_beta,
              p_vals = sigset_merge$final_p_val)

  return(out)
}



