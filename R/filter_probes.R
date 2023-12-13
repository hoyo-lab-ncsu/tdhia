

#' filter_probes
#'
#' Applies QC and filtering for a probe_beta matrix.
#'
#' @param probe_beta a list object that a dataframe of beta values and pvalues.
#'   $probe_beta_df: a dataframe of beta values, probe_id x sample_id
#'   $probe_pval_df: a dataframe of signal detection p-values, probe_id x sample_id
#' @param discard_unmapped_probes boolean flag, when TRUE probes that do not map
#' to a unique position in the genome are discarded from the probe beta and
#' p-value dataframes.
#' @param max_sig_pval a maximum threshold for the signal detection value, beta
#' values higher than this max are set to NA.
#' @param max_probe_fail_rate proportion value [0-1], if a probe fails the p-value
#' threshold above this proportion of the samples, all measurements from this probe
#' are set to NA.
#' @param discard_failed_probes boolean flag, when TRUE any probes that have NA
#' values across all samples are discarded.
#' @param mft manifest file for particular imprintome array used in data collection.
#' A dataframe with metadata mapping probe_ids to CpG sites. See ?manifest_v1A2
#' for more information on this metadata.
#'
#' @returns a dataframe containing beta values where rows are probe_ids and columns
#' are the sample_ids (basenames of the IDAT files or some other mapping
#' specified by idat_remappings).
#'
filter_probes <- function(probe_beta, discard_unmapped_probes = TRUE,
                          max_sig_pval = 0.2, max_probe_fail_rate = 0.5,
                          discard_failed_probes = TRUE) {
  # save(list = ls(all.names = TRUE), file = "filter_probes_debug.RData")
  # load(file = "filter_probes_debug.RData")

  # Load manifest
  mft = probe_beta$manifest


  # Discard probe_ids(s) that do not map uniquely to a single genome location
  #----------------------------------------------------------------------------
  #   (Defined in manifest file where MAPINFO == (0 or NA))
  #   Data is copied over to filter_* probe_beta and probe_ pval matrices.
  if (discard_unmapped_probes) {
    f_unmmaped = function (x) is.na(x) | x==0;
    unmapped_probe_ids <- sapply(mft$Probe_ID,
                                 function(x) f_unmmaped(mft$MAPINFO[which(x==mft$Probe_ID)[1]]))
    filt_probe_beta_df <- probe_beta$probe_beta_df[!unmapped_probe_ids,]
    filt_probe_pval_df <- probe_beta$probe_pval_df[!unmapped_probe_ids,]
  } else {
    # Copy the probe_beta and p-value matrices and beginning filtering steps
    filt_probe_beta_df <- probe_beta$probe_beta_df
    filt_probe_pval_df <- probe_beta$probe_pval_df
  }


  # Set individual probe_beta values to NA if the p-value is above max threshold
  #-----------------------------------------------------------------------------
  if (!is.null(max_sig_pval)) {
    filt_probe_beta_df[filt_probe_pval_df < max_sig_pval] <- NA
  }


  # Set entire row of probe measures to NA if not enough pass p-value check above
  #-----------------------------------------------------------------------------
  if(!is.null(max_probe_fail_rate)) {
    # Binary index of rows, TRUE means row failed QC
    row_fail_ind <- rowSums(is.na(filt_probe_beta_df))/ncol(filt_probe_beta_df) <
      max_probe_fail_rate
    # Set rows to NA
    filt_probe_beta_df[row_fail_ind,] <- NA
  }


  # Discard probes/rows where all samples have a NA for the beta value
  #-----------------------------------------------------------------------------
  if(discard_failed_probes) {
    # Binary index of rows, TRUE means row failed QC
    row_all_na_ind <- rowSums(is.na(filt_probe_beta_df)) == ncol(filt_probe_beta_df)
    # Set rows to NA
    filt_probe_beta_df <- filt_probe_beta_df[!row_fail_ind,]
    filt_probe_pval_df <- filt_probe_pval_df[!row_fail_ind,]
  }


 probe_beta <- list(probe_beta_df = filt_probe_beta_df,
                      probe_pval_df = filt_probe_pval_df,
                      platform = probe_beta$platform,
                      manifest = probe_beta$manifest
  )

  return(probe_beta)


}



