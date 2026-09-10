

#' tdhia_pipeline
#'
#' @description runs standard tdhia pipeline in a single step.
#'
#' TODO: this function does not check if the cached data with probe_data matches
#'  the settings of the current function cell, need to add in future.
#' @param idat_dir_paths a vector of strings specifying full directory paths
#' where IDAT files are located.
#' @param probe_data_cache either a probe_data object, a filepath to processed data, or NULL.
#'   1) probe_data: the initial processing of idat files to probe_data is
#'    skipped. All donstream processing still occurs.
#'   2) filepath: a path where probe data is saved. If the file exists, the 
#'   probe data is loaded from there instead of reprocessing.
#'   3) NULL: the probe data is processed and not cached to disk.
#' @param multicore boolean, when true uses n.cores processing in sesame.
#' @param idat_basenames vector of basenames for idat files, may be a subset of 
#' files that are stored in the cache saved in probe_data_cache.
#' @param discard_unmapped_probes boolean, when true, probe_ids
#' @param max_sig_pval single numeric between (0-1), maximum sesame signal
#' p-value that is allowed
#' @param set_failed_betas_na boolean, when true probe measurements that fail
#' the p-value threshold are set to NA (missing).
#' @param max_probe_fail_rate maximum number of measurements that a probe can
#' fail the p-value threshold before being removed entirely from dataset (all
#' values set to NA).
#' @param discard_failed_probes boolean when true, the probes that fail the
#' max_probe_fail_rate are removed entirely.
#' @param db_flag flag for saving function work space for debugging.
#' @param smooth_adj_cpgs boolean, when true, adjacent cpg beta values within an
#' ICR are smoothed with a sliding window average.
#' @param max_icr_fail_rate maximum allowable fraction of samples that can be
#' missing before the entire ICR is excluded from dataset.
#' @param overwrite_probe_cache todo
#' @param max_patient_fail_rate maximum number of imprintome cpg sites that fail
#' the max_sig_pval threshold. Any patients with a higher fraction will be 
#' discarded.
#' @param merge_replicates string with the following possible values that
#' determines how replicate probes are handled.
#'  - NULL: the probe replicates are not merged (default value).
#'  - pre_beta: replicates are merged by averaging fluorescent signal from
#'   each channel individually, beta values are then calculated from these
#'   averaged values (pre_beta = merging done before beta calculation).
#'  - post_beta: beta value is calculated before replicates are merged and beta
#'   values are averaged between replicates (post_beta = merging done before
#'   beta calculation).
#' @param enforce_req_idats check that all idat files requested in idat_basenames
#'  argument are found on disk. Throws error if this is not the case.
#' @param min_design_score discard probes below this design score threshold.
#' Default: NA (don't discard any probes by design score). Design scores vary
#' from 0-1, with a higher value being a better designed probe.
#' @param probe_beta provide probe data from a previously processed function call
#'  of load_idata_to_probes. This is useful to provide when you want different
#'   settings downstead of probe beta culations (such as setting individual beta 
#'   values to NA if the sesame p-valu threshold is not met).
#' @return named list with all of the output from each of the pipeline steps.
#' @export
tdhia_pipeline <- function(
    idat_dir_paths = NULL, overwrite_probe_cache = F, probe_data_cache = NULL,
    multicore = TRUE, idat_basenames = NULL, discard_unmapped_probes = TRUE , 
    max_sig_pval = 0.2, set_failed_betas_na = FALSE, max_probe_fail_rate = 0.2,
    discard_failed_probes = TRUE, max_patient_fail_rate = 0.25, smooth_adj_cpgs = FALSE,
    max_icr_fail_rate = 0.2, db_flag = FALSE, merge_replicates = "pre_beta",
    enforce_req_idats = TRUE, min_design_score = NA, verbose = T) {

  if(db_flag) save(list = ls(all.names = TRUE), file = "tdhia_pipeline.RData")
  # load(file = "tdhia_pipeline.RData")
  if (is.character(probe_data_cache)) dir.create(dirname(probe_data_cache),
                                                 showWarnings = F, recursive = T)
  
  verbosecat = \(x) if (verbose) cat(x)
  
  # Store all data in fields of list
  data_beta = list()
  
  # 1) process idat files to probe data
  if ( is.null(probe_data_cache) ||
       (is.character(probe_data_cache) && !file.exists(probe_data_cache)) ||
       (is.list(probe_data_cache) && overwrite_probe_cache) ) {
    # If probe_data_cache is null, OR is a path and does not exist, then process
    verbosecat("Processing probe data.\n")
    data_beta$probe_beta <-
      load_idata_to_probes(idat_dir_paths = idat_dir_paths, multicore = multicore,
                           idat_basenames = idat_basenames, 
                           quantile_norm = FALSE, db_flag = db_flag, 
                           merge_replicates = merge_replicates, 
                           enforce_req_idats = enforce_req_idats)
    
    # Save data to disk if probe_data_cache is a path
    if (is.character(probe_data_cache)) saveRDS(data_beta$probe_beta, file = probe_data_cache)
    
  } else if (is.list(probe_data_cache)) {
    # If probe_data_cache is a probe_data object, then assign
    verbosecat("Using processed probe data from input argument.\n")
    data_beta$probe_beta <- probe_data_cache
  } else if (is.character(probe_data_cache)) { 
    # Load if probe_data_cache is a path
    verbosecat("Loading processed probe data from cache file.\n")
    data_beta$probe_beta <- readRDS(probe_data_cache)
  } else { stop("probe_beta_cache logic flawed, something went wrong")}

  # Ensure that requested idats match the ones in probe beta
  if (!is.null(idat_basenames)) {
    data_beta$probe_beta <- subset_probe_beta(data_beta$probe_beta, idat_basenames)
  }
  
  # 2) Filter probes that are not mapped and discard poor signal
  verbosecat("Filtering probe data.\n")
  data_beta$filt_probe_beta <- filter_probes(
    probe_beta = data_beta$probe_beta, discard_unmapped_probes = discard_unmapped_probes,
    max_sig_pval = max_sig_pval, set_failed_betas_na = set_failed_betas_na, 
    max_probe_fail_rate = max_probe_fail_rate, min_design_score = min_design_score,
    discard_failed_probes = discard_failed_probes, max_patient_fail_rate = max_patient_fail_rate, 
    db_flag = db_flag)
  
  
  #3  Convert probe beta matrix to a cpg beta matrix
  verbosecat("Processing probe data to cpg data\n")
  data_beta$cpg_beta <- convert_probes_to_cpgs(
    data_beta$filt_probe_beta, quantile_norm = FALSE,  db_flag = db_flag, 
    smooth_adj_cpgs = smooth_adj_cpgs)
  
  
  #4  Convert probe beta matrix to an icr beta matrix
  verbosecat("Processing cpg data to icr data\n")
  data_beta$icr_beta <- convert_cpgs_to_icrs(data_beta$cpg_beta, 
                                             max_icr_fail_rate = max_icr_fail_rate)
  
  # document input args for reproducibility
  data_beta$input_args = input_args = mget(names(formals()), envir = environment())
  
  return(data_beta)
  
}



subset_probe_beta <- function(probe_beta, idat_basenames) {
  # Update the beta, p-value matrices and list of basenames
  probe_beta$probe_beta_df <-
    probe_beta$probe_beta_df %>% select(all_of(idat_basenames))
  probe_beta$probe_pval_df <-
    probe_beta$probe_pval_df %>% select(all_of(idat_basenames))
  probe_beta$idat_filepaths <- intersect(probe_beta$idat_filepaths,  idat_basenames)
  probe_beta$input_args$idat_basenames <-idat_basenames
  
  return(probe_beta)
  
}