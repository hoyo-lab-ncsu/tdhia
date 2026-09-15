

#' Run the Imprintome Processing Pipeline
#'
#' Loads IDAT measurements or cached probe data, filters probes and samples,
#' and aggregates methylation beta values to CpG and ICR levels.
#'
#' @param idat_dir_paths Character vector of directories containing IDAT files,
#'   passed to load_idata_to_probes() when probe data must be processed.
#' @param overwrite_probe_cache Logical; reprocess an input probe-data list
#'   when TRUE. An existing cache file supplied as a path is still reused.
#' @param probe_data_cache NULL, a probe-data list from load_idata_to_probes(),
#'   or an RDS path. NULL processes without saving a cache; a missing path is
#'   populated and an existing path is read. Parent directories are created.
#' @param multicore Logical or core count passed to load_idata_to_probes().
#' @param idat_basenames Character vector of requested sample basenames.
#'   NULL uses all available samples; supplied names also subset cached data.
#' @param discard_unmapped_probes Logical; retain only mapped CpG probes when
#'   TRUE, as determined by filter_probes().
#' @param max_sig_pval Numeric detection p-value threshold. Measurements pass
#'   when their detection p-value is strictly below this value.
#' @param set_failed_betas_na Logical; replace individual failed measurements
#'   with NA. Probe-level and sample-level filtering still run when FALSE.
#' @param max_probe_fail_rate Maximum fraction of failed measurements per
#'   probe. Probes above the threshold have all beta values set to NA.
#' @param discard_failed_probes Logical; remove probes whose beta values are
#'   all missing after filtering.
#' @param max_patient_fail_rate Maximum failure fraction per sample, passed
#'   to filter_probes(); failing samples have all beta values set to NA and
#'   are removed using that helper's default discard_failed_patients setting.
#' @param smooth_adj_cpgs Logical; apply a centered three-CpG rolling mean
#'   within each ICR during CpG aggregation.
#' @param max_icr_fail_rate Numeric threshold passed to convert_cpgs_to_icrs().
#'   Its current implementation reports failures but does not apply that
#'   threshold to the returned ICR table.
#' @param db_flag Logical; save debugging workspaces for this function and
#'   processing helpers that receive the flag.
#' @param merge_replicates Replicate handling passed to the IDAT loader:
#'   "pre_beta" (default) averages fluorescence before calculating beta,
#'   "post_beta" averages beta values, and NULL leaves replicates separate.
#' @param enforce_req_idats Logical; require all requested IDAT basenames
#'   when loading files.
#' @param min_design_score Minimum probe design score. NA disables this
#'   filter; otherwise only scores strictly greater than the threshold remain.
#' @param verbose Logical; print pipeline progress. This flag is not forwarded
#'   to the processing helpers, which can print their own messages.
#'
#' @details
#' Filtering and both aggregation steps run even when probe data are cached.
#' Cached probe data are not checked against the current processing settings.
#' Quantile normalization is disabled in the loader and CpG conversion calls.
#' See the individual processing functions for filtering and aggregation
#' details.
#'
#' @return A named list with:
#'   - probe_beta: loaded or cached probe data, subset to requested samples.
#'   - filt_probe_beta: filtered probe data and QC summaries from
#'     filter_probes().
#'   - cpg_beta: CpG beta values and metadata from convert_probes_to_cpgs().
#'   - icr_beta: ICR beta values and metadata from convert_cpgs_to_icrs().
#'   - input_args: argument values captured for reproducing the pipeline call.
#' @seealso [load_idata_to_probes()], [filter_probes()],
#'   [convert_probes_to_cpgs()], [convert_cpgs_to_icrs()]
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