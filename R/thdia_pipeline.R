

#' tdhia_pipeline
#'
#' @description runs standard tdhia pipeline in a single step.
#'
#' @param idat_dir_paths a vector of strings specifying full directory paths
#' where IDAT files are located.
#' @param multicore boolean, when true uses n.cores processing in sesame.
#' @param idat_basenames rename idat_basenames to a patient/ study id
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
#'
#' @return named list with all of the output from each of the pipeline steps.
#'
tdhia_pipeline <- function(idat_dir_paths = NULL, multicore = TRUE,
                           idat_basenames = NULL,
                           discard_unmapped_probes = TRUE , max_sig_pval = 0.2,
                           set_failed_betas_na = FALSE, max_probe_fail_rate = 0.25,
                           discard_failed_probes = TRUE, smooth_adj_cpgs = FALSE,
                           max_icr_fail_rate = 0.2, db_flag = FALSE) {


  if(db_flag) save(list = ls(all.names = TRUE), file = "tdhia_pipeline.RData")
  # load(file = "tdhia_pipeline.RData")

  # 1) Load IDATS and convert to probe beta matrix
  # Output: probe_id(row) x patient(col)
    probe_beta <-
      load_idata_to_probes(idat_dir_paths = idat_dir_paths, multicore = multicore,
                           idat_basenames = idat_basenames,
                           quantile_norm = FALSE)

  # 2) Filter probes that are not mapped and discard poor signal
  # Output: probe_id(row) x patient(col)
  filt_probe_beta <-
    filter_probes(probe_beta = probe_beta, discard_unmapped_probes = discard_unmapped_probes,
                  max_sig_pval = max_sig_pval, set_failed_betas_na = set_failed_betas_na,
                  max_probe_fail_rate = max_probe_fail_rate,
                  discard_failed_probes = discard_failed_probes, db_flag = db_flag)

  #3  Convert probe beta matrix to a cpg beta matrix
  # Output: cpg_id(row) x patient(cols)
  cpg_beta <- convert_probes_to_cpgs(filt_probe_beta, quantile_norm = FALSE,
                                     db_flag = db_flag, smooth_adj_cpgs = smooth_adj_cpgs)

  #4  Convert probe beta matrix to an icr beta matrix
  # Output: icr_id(row) x patient(cols)
  icr_beta <- convert_cpgs_to_icrs(cpg_beta, max_icr_fail_rate = max_icr_fail_rate)

  # Assemble Output structure
  out <- list(probe_beta = probe_beta, filt_probe_beta = filt_probe_beta,
              cpg_beta = cpg_beta,
       icr_beta = icr_beta)
  return(out)
}
