

#' filter_probes
#'
#' Applies probe and measurement filtering for a probe_beta matrix.
#'
#' @param probe_beta a list object with the following named fields
#'   $probe_beta_df: a dataframe of beta values, probe_id x sample_id
#'   $probe_pval_df: a dataframe of signal detection p-values, probe_id x sample_id
#' @param discard_unmapped_probes boolean flag, when TRUE probes that do not map
#' to a unique genomic position or CpG site are discarded from the probe beta and
#' p-value dataframes. Probes that match to a CpG site are identified by their prefix
#' in the probe name (<prefix><index>_<version tag>). Prefixes are (ctl:control,
#' cg:CpG, ch:CHG, mu:multi-unique, rp:repetitive element, or rs:SNP)
#' @param max_sig_pval a maximum threshold for the signal detection value, beta
#' values higher than this max are set to NA.
#' @param set_failed_betas_na boolean flag, when TRUE beta measurements with pval >
#' max_sig_pval will be set to NA. When false, the beta value is left as is.
#' @param max_probe_fail_rate proportion value (0-1), if a probe fails the p-value
#' threshold above this proportion of the samples, all measurements from this probe
#' are set to NA.
#' @param discard_failed_probes boolean flag, when TRUE any probes that have NA
#' values across all samples are discarded (default = TRUE)
#' @param min_design_score discard probes below this design score threshold.
#' Default: NA (don't discard any probes by design score). Design scores vary
#' from 0-1, with a higher value being a better designed probe.
#' @param icr_conf_levels a numerical vector listing the desired ICR confidence 
#' levels for filtered probes.Any probes that are not included in this list are 
#' excluded. ICR confidence scores are (1-3), where
#' (1) High: literature validated gold standard.
#' (2) Medium: experimental evidence of gametic origin of differential methylation.
#' (3) Low: experimental evidence of 50% methylation in somatic tissues.
#' Example: icr_conf_levels = c(1,2)
#'   Would discard any probes associated with a ICR confidence score of 3.
#' @param verbose boolean flag, when TRUE prints the results of each filtering
#' step (default = TRUE).
#' @param db_flag boolean when true export workspace to disk for debugging.
#' @param mft manifest dataframe for the array if nto specified, the default is
#' loaded.
#'
#' @returns a dataframe containing beta values where rows are probe_ids and columns
#' are the sample_ids (basenames of the IDAT files or some other mapping
#' specified by idat_remappings).
#'
#' @importFrom rlang .data
#' @export
filter_probes <- function(probe_beta, discard_unmapped_probes = TRUE,
                          max_sig_pval = 0.2, set_failed_betas_na = TRUE, mft = NULL,
                          max_probe_fail_rate = 0.5, discard_failed_probes = TRUE,
                          min_design_score = NA, icr_conf_levels = NULL,
                          verbose = TRUE, db_flag = FALSE) {
  if (db_flag) {save(list = ls(all.names = TRUE), file = "filter_probes_debug.RData")}
  # load(file = "filter_probes_debug.RData")

  # Only print cat() output if the user requests it
  verbosecat <-function(x) if (verbose) cat(x)

  # Extract manifest dataframe
  if (is.null(mft)) {
  cat("Loading manifest from internal data\n"); mft = tdhia::manifest_v1A2_design_scores
  }
  verbosecat(sprintf("Probe manifest: manifest file has a total of %.0f probes.\n\n",
              nrow(mft)))

  # Report probes removed from previous steps                               ####
  #_____________________________________________________________________________
  verbosecat(sprintf("Probe filter: %.0f%% probes discarded ( %i/ %i) from previous steps
                %.0f probes now remain.\n\n",
                     100*(nrow(mft) - nrow(probe_beta$probe_beta_df))/nrow(mft),
                     nrow(mft)- nrow(probe_beta$probe_beta_df),
                     nrow(mft),
                     nrow(probe_beta$probe_beta_df)))

  # Discard probe_ids(s) that do not map uniquely to genomic location and CpG ####
  #_____________________________________________________________________________
  #   (Defined in manifest file where MAPINFO == (0 or NA))
  #   Data is copied over to filter_* probe_beta and probe_ pval matrices.
  if (discard_unmapped_probes) {
    # Returns true of probe is no genomic location
    f_unmmaped = function (x) is.na(x) | x==0;
    f_match = function (x, probe_ids) { match_ind = which(x == mft$Probe_ID)
      if(length(match_ind)==0) match_ind = which(stringr::str_replace(x, "_.{4}$","") ==
                                                   mft$cpg_id)[1]
      return(match_ind)
    }
    # unmapped_probe_ids <- sapply(mft$Probe_ID,
    #                              function(x) f_unmmaped(mft$MAPINFO[f_match(x,mft$Probe_ID)[1]]))


    # unmapped_probe_ids
    mapped_probe_ids <- sapply(rownames(probe_beta$probe_beta_df),
                                 function(x) !f_unmmaped(mft$MAPINFO[f_match(x,mft$Probe_ID)[1]]))

    #
    # probe_beta %>% select()
    # left_join(x = data.frame(cpg_id = rownames(probe_beta$probe_beta_df)),
    #           y = 1)


    # Find probes that target a cg site
    cpg_ids <- grepl("^cg[0-9]",rownames(probe_beta$probe_beta_df))
    # For reference probe_id stars with either:
    # ctl: control, cg: CPG, ch: CHG, mu: multi-unique,
    # rp: repetitive element, or rs: SNP probe

    # Keep probes that are mapped to a cpg site in the imprintome
    filt_probe_beta_df <- probe_beta$probe_beta_df[mapped_probe_ids & cpg_ids,]
    filt_probe_pval_df <- probe_beta$probe_pval_df[mapped_probe_ids & cpg_ids,]

    verbosecat(sprintf("Probe filter: discarding %.0f%% probes ( %i/ %i) in dataset b/c they don't
                map uniquely to the genome or a CpG site.
                %.0f probes now remain.\n\n",
                100*nrow(filt_probe_beta_df)/nrow(probe_beta$probe_beta_df),
                nrow(probe_beta$probe_beta_df)- nrow(filt_probe_beta_df),
                nrow(probe_beta$probe_beta_df),
                nrow(filt_probe_beta_df)))

  } else {
    # Copy the probe_beta and p-value matrices and beginning filtering steps
    filt_probe_beta_df <- probe_beta$probe_beta_df
    filt_probe_pval_df <- probe_beta$probe_pval_df
  }


  # Set individual probe_beta values to NA if the p-value < threshold ####
  #_____________________________________________________________________________
  sig_pval_pass <- filt_probe_pval_df < max_sig_pval
  if ((!is.null(max_sig_pval)) && set_failed_betas_na) {
    filt_probe_beta_df[!sig_pval_pass] <- NA
    verbosecat(sprintf("Probe filter: %.0f%% of probe beta measurements ( %i/ %i) failed
                 the signal max p-value threshold of %.2f, setting them to NA.\n\n",
                100*sum(is.na(filt_probe_beta_df))/prod(dim(filt_probe_beta_df)),
                sum(is.na(filt_probe_beta_df)),prod(dim(filt_probe_beta_df)),
                max_sig_pval))

    # Five number summary of failure rates
    verbosecat('Distribution of failiure rates for probes:\n')
    qsum <- (stats::quantile(unname(rowSums(is.na(filt_probe_beta_df))/ncol(filt_probe_beta_df)),
             seq(.1,1,.1)))
    df_sum = data.frame(Quantile = names(qsum), "Fraction.Failed" = unname(qsum))
    # print(df_sum)
    for (n in 1:nrow(df_sum)) {
      verbosecat(sprintf("%s \t %.03f\n", df_sum$Quantile[n], df_sum$Fraction.Failed[n]))
    }
    verbosecat("\n")
    
  } else {
    verbosecat(sprintf("Probe filter: keeping all individual beta values that
                       exceed p-value.\n\n"))
  }
  # hist(rowSums(!is.na(filt_probe_beta_df))/ncol(filt_probe_beta_df),
  #      main = paste("Histogram of probe pass rate"),
  #      xlab = "Fraction passed probes")


  # Set entire row of probe measures to NA if not enough pass p-value check ####
  #_____________________________________________________________________________
  if(!is.null(max_probe_fail_rate)) {
    # Binary index of rows, TRUE means row failed QC
    probe_fail_rate_df <- data.frame(cpg_id = rownames(filt_probe_beta_df),
                                probe_fail_rate = rowSums(!sig_pval_pass) /
                                  ncol(filt_probe_beta_df))

    row_fail_flag <-  probe_fail_rate_df$probe_fail_rate > max_probe_fail_rate
    # Set rows to NA
    filt_probe_beta_df[row_fail_flag,] <- NA
    verbosecat(sprintf("Probe filter: %.0f%% of probes ( %i/ %i) had a signal p-value
                fail rate above threshold of %.0f%%, setting all beta values to NA
                for those probes.\n\n",
        100*sum(row_fail_flag)/length(row_fail_flag), sum(row_fail_flag),
        length(row_fail_flag), 100*max_probe_fail_rate))
  } else {probe_fail_rate_df = NULL}
  # hist(rowSums(!is.na(filt_probe_beta_df))/ncol(filt_probe_beta_df),
  #      main = paste("Histogram of probe pass rate"),
  #      xlab = "Fraction passed probes")


  # Discard probes/rows where all samples have a NA for the beta value   ######
  #_____________________________________________________________________________
  if(discard_failed_probes) {
    # Binary index of rows, TRUE means row failed QC
    na_keep_ind <- rowSums(is.na(filt_probe_beta_df)) != ncol(filt_probe_beta_df)
    # Set rows to NA
    filt_probe_beta_df <- filt_probe_beta_df[na_keep_ind, ]
    filt_probe_pval_df <- filt_probe_pval_df[na_keep_ind, ]
    verbosecat(sprintf("Probe filter: discarding %.0f%% probes ( %i/ %i) because all measurements are now NA.
              %i probes now remain. \n\n",
                       100*sum(na_keep_ind==FALSE)/length(na_keep_ind),
                       sum(na_keep_ind==FALSE),
                       length(na_keep_ind),
                       sum(na_keep_ind==TRUE)))
  }



  # Filter probes with low design scores              ##########################
  #_____________________________________________________________________________
  design_score_df <- dplyr::left_join(
    x = data.frame(Name = stringr::str_replace(rownames(filt_probe_beta_df),"_.{4}$","")),
    y = dplyr::select(mft, c("Name", "Design.Score")), by = "Name",
    unmatched = "drop", keep = FALSE, multiple = "first")
  design_score_keep_flag = design_score_df$Design.Score > min_design_score

  if (!is.na(min_design_score)) {
    filt_probe_beta_df <- filt_probe_beta_df[design_score_keep_flag, ]
    filt_probe_pval_df <- filt_probe_pval_df[design_score_keep_flag, ]
    design_score_df <- design_score_df[design_score_keep_flag, ]

    verbosecat(sprintf("Probe filter: discarding %.2f%% probes ( %i/ %i) because of low design score.
              %i probes now remain. \n\n",
                       100*sum(design_score_keep_flag==FALSE)/length(design_score_keep_flag),
                       sum(design_score_keep_flag==FALSE),
                       length(design_score_keep_flag),
                       sum(design_score_keep_flag==TRUE)))
  }

   
  # Filter probes by ICR confidence level                             ##########
  #_____________________________________________________________________________
  if (!is.null(icr_conf_levels)) {
    # icr_conf_levels is considered a WHITELIST, not a threshold.
    df_metadata <- add_metadata_from_cpg(str_replace(rownames(filt_probe_beta_df), "_.*$",""))
    icr_conf_keep_flag <- df_metadata$icr_conf %in% icr_conf_levels
    
    # FIlter out probes that are not part of specified icr_conf levels
    filt_probe_beta_df <- filt_probe_beta_df[icr_conf_keep_flag, ]
    filt_probe_pval_df <- filt_probe_pval_df[icr_conf_keep_flag, ]
    design_score_df <- design_score_df[icr_conf_keep_flag, ]
    
    verbosecat(sprintf("Probe filter: discarding %.2f%% of probes ( %i/ %i) from undesired ICR confidence levels.
              %i probes now remain. \n\n",
                       100*sum(icr_conf_keep_flag==FALSE)/length(icr_conf_keep_flag),
                       sum(icr_conf_keep_flag==FALSE),
                       length(icr_conf_keep_flag),
                       sum(icr_conf_keep_flag==TRUE)))
  }
  
  
  # Report how many missing measurements at end of all filtering steps.
  verbosecat(sprintf("Probe filter: After all filtering, %.0f%% of individual probe measurments ( %i/ %i) are now NA.\n\n",
                     100*sum(is.na(filt_probe_beta_df)) / prod(dim(filt_probe_beta_df)),
                     sum(is.na(filt_probe_beta_df)),
                     prod(dim(filt_probe_beta_df))))

  # Report how many missing measurements at end of all filtering steps.
  verbosecat(sprintf("Probe filter: After all filtering, %.0f%% of probes ( %i/ %i) remain.\n\n",
                     100*nrow(filt_probe_beta_df) / nrow(probe_beta$probe_beta_df),
                     nrow(filt_probe_beta_df),
                     nrow(probe_beta$probe_beta_df)))

  probe_beta <- list(probe_beta_df = filt_probe_beta_df,
                     probe_pval_df = filt_probe_pval_df,
                     design_scores = design_score_df,
                     platform = probe_beta$platform,
                     manifest = probe_beta$manifest,
                     probe_fail_rate_df = probe_fail_rate_df
  )

  return(probe_beta)


}




#' pbates
#'
#'  Cumulative distribution function for the Bates distribution, which is the mean
#' of n observations sampled from a uniform distribution between \[0,1]. Uses a
#' simple Monte Carlo approach to model the null distribution. This can be used
#' to calculate a new signal detection p-value across replicate probes that are
#' found in the imprintome array. The random variable would be the mean p-value
#' from replicate probes.This function is most efficient if all mean-p-values are
#' processed at once because it pre-catches the ECDF from a bates distribution
#' with different sample sizes.
#'
#' @param mean_p_vals vector of observed means from a series of samples
#' @param samples vector of the number of observations per sample
#' @param n_sims number of simulaations
#'
#' @export
pbates <- function(mean_p_vals, samples, n_sims = 1e6) {

  # Vectorize(pbates_unvector,  vectorize.args = c("obs_pval_bar", "n"))
  p_vals = rep(NA, length(mean_p_vals))


  # Pre-generate ECDFs for bates dstribution with different sample counts
  unq_n <- unique(samples)
  unq_n <- unq_n[unq_n>1]
  ecdfs <- list()
  for (k in seq_along(unq_n)) {
    ecdfs[[as.character(unq_n[k])]] = stats::ecdf(
      Matrix::rowSums(matrix(dqrng::dqrunif(unq_n[k]*n_sims, 0, 1),
                             ncol = unq_n[k]))/unq_n[k])
  }

  # Calculate bates p-values using the ECDF that has the same sample count
  # For cases with only 1 sample, then return original p-value
  for (k in seq_along(mean_p_vals)) {
    if (samples[k]==1) {
      p_vals[k] = mean_p_vals[k]
    } else {
      p_vals[k] = ecdfs[[as.character(samples[k])]](mean_p_vals[k])
    }
  }

  return(p_vals)
}


