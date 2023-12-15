

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
#' @param max_probe_fail_rate proportion value (0-1), if a probe fails the p-value
#' threshold above this proportion of the samples, all measurements from this probe
#' are set to NA.
#' @param discard_failed_probes boolean flag, when TRUE any probes that have NA
#' values across all samples are discarded.
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

  # Extract manifest dataframe
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

    cat(sprintf("Probe filtering: discarding %.0f%% probes ( %i/ %i) b/c they don't
                map to genome. %.0f probes now remain.\n",
                100*nrow(filt_probe_beta_df)/nrow(probe_beta$probe_beta_df),
                nrow(probe_beta$probe_beta_df)- nrow(filt_probe_beta_df),
                nrow(probe_beta$probe_beta_df),
                nrow(filt_probe_beta_df)))
  } else {
    # Copy the probe_beta and p-value matrices and beginning filtering steps
    filt_probe_beta_df <- probe_beta$probe_beta_df
    filt_probe_pval_df <- probe_beta$probe_pval_df
  }


  # Set individual probe_beta values to NA if the p-value is above max threshold
  #-----------------------------------------------------------------------------
  if (!is.null(max_sig_pval)) {
    filt_probe_beta_df[filt_probe_pval_df > max_sig_pval] <- NA
    cat(sprintf("Probe filtering: %.0f%% of beta measurements ( %i/ %i) failed the
                signal p-value max threshold of %.2f, setting them to NA.\n",
        100*sum(is.na(filt_probe_beta_df))/prod(dim(filt_probe_beta_df)),
        sum(is.na(filt_probe_beta_df)),prod(dim(filt_probe_beta_df)),
        max_sig_pval))
  }
  # hist(rowSums(!is.na(filt_probe_beta_df))/ncol(filt_probe_beta_df),
  #      main = paste("Histogram of probe pass rate"),
  #      xlab = "Fraction passed probes")


  # Set entire row of probe measures to NA if not enough pass p-value check above
  #-----------------------------------------------------------------------------
  if(!is.null(max_probe_fail_rate)) {
    # Binary index of rows, TRUE means row failed QC
    row_fail_flag <- rowSums(is.na(filt_probe_beta_df))/ncol(filt_probe_beta_df) >
      max_probe_fail_rate
    # Set rows to NA
    filt_probe_beta_df[row_fail_flag,] <- NA
    cat(sprintf("Probe filtering: %.0f%% of probes ( %i/ %i) had a signal p-value
                fail rate above %.0f%%, setting all beta values to NA for those probes.\n",
        100*sum(row_fail_flag)/length(row_fail_flag), sum(row_fail_flag),
        length(row_fail_flag), 100*max_probe_fail_rate))
  }
  # hist(rowSums(!is.na(filt_probe_beta_df))/ncol(filt_probe_beta_df),
  #      main = paste("Histogram of probe pass rate"),
  #      xlab = "Fraction passed probes")


  # Discard probes/rows where all samples have a NA for the beta value
  #-----------------------------------------------------------------------------
  if(discard_failed_probes) {
    # Binary index of rows, TRUE means row failed QC
    row_all_na_ind <- rowSums(is.na(filt_probe_beta_df)) == ncol(filt_probe_beta_df)
    # Set rows to NA
    filt_probe_beta_df2 <- filt_probe_beta_df[!row_all_na_ind,]
    filt_probe_pval_df2 <- filt_probe_pval_df[!row_all_na_ind,]
  } else {
    filt_probe_beta_df2 <- filt_probe_beta_df
    filt_probe_pval_df2 <- filt_probe_pval_df
  }
  cat(sprintf("Probe filtering: discarded %.0f%% probes ( %i/ %i) because all measurements were NA.
              %i probes now remain. \n",
      100*(nrow(filt_probe_beta_df) - nrow(filt_probe_beta_df2))/nrow(filt_probe_beta_df),
      nrow(filt_probe_beta_df) - nrow(filt_probe_beta_df2),
      nrow(filt_probe_beta_df), nrow(filt_probe_beta_df2)))

  # hist(rowSums(!is.na(filt_probe_beta_df2))/ncol(filt_probe_beta_df2),
  #      main = paste("Histogram of probe pass rate"),
  #      xlab = "Fraction passed probes")


 probe_beta <- list(probe_beta_df = filt_probe_beta_df,
                      probe_pval_df = filt_probe_pval_df,
                      platform = probe_beta$platform,
                      manifest = probe_beta$manifest
  )

  return(probe_beta)


}



