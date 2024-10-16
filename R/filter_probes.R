

#' filter_probes
#'
#' Applies QC and filtering for a probe_beta matrix.
#'
#' @param probe_beta a list object that a dataframe of beta values and pvalues.
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
#' @param min_probe_score discard probes below this design score threshold.
#' @param verbose boolean flag, when TRUE prints the results of each filtering
#' step (default = TRUE).
#' @param db_flag boolean when true export workspace to disk for debugging.
#'
#' @returns a dataframe containing beta values where rows are probe_ids and columns
#' are the sample_ids (basenames of the IDAT files or some other mapping
#' specified by idat_remappings).
#'
#' @importFrom rlang .data
#' @export
filter_probes <- function(probe_beta, discard_unmapped_probes = TRUE,
                          max_sig_pval = 0.2, set_failed_betas_na = TRUE,
                          max_probe_fail_rate = 0.5, min_probe_score = 0.2,
                          discard_failed_probes = TRUE, verbose = TRUE,
                          db_flag = FALSE) {
  if (db_flag) {save(list = ls(all.names = TRUE), file = "filter_probes_debug.RData")}
  # load(file = "filter_probes_debug.RData")

  # Only print cat() output if the user requests it
  verbosecat <-function(x) if (verbose) cat(x)

  # Extract manifest dataframe
  if (!is.null(probe_beta$manifest)) {
  cat("Loading manifest from probe_beta data\n");  mft = probe_beta$manifest
  } else {
  cat("Loading manifest from internal data\n"); mft = tdhia::manifest_v1A2
  }
  verbosecat(sprintf("Probe manifest: manifest file has a total of %.0f probes.\n\n",
              nrow(mft)))


  # Discard probe_ids(s) that do not map uniquely to genomic location and CpG
  #----------------------------------------------------------------------------
  #   (Defined in manifest file where MAPINFO == (0 or NA))
  #   Data is copied over to filter_* probe_beta and probe_ pval matrices.
  if (discard_unmapped_probes) {
    f_unmmaped = function (x) is.na(x) | x==0;
    unmapped_probe_ids <- sapply(mft$Probe_ID,
                                 function(x) f_unmmaped(mft$MAPINFO[which(x==mft$Probe_ID)[1]]))

    # Find probes that target a cg site
    cpg_ids <- grepl("^cg",rownames(probe_beta$probe_beta_df))
    # For reference probe_id stars with either:
    # ctl: control, cg: CPG, ch: CHG, mu: multi-unique,
    # rp: repetitive element, or rs: SNP probe

    # Keep probes that are mapped to a cpg site in the imprintome
    filt_probe_beta_df <- probe_beta$probe_beta_df[!unmapped_probe_ids & cpg_ids,]
    filt_probe_pval_df <- probe_beta$probe_pval_df[!unmapped_probe_ids & cpg_ids,]

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


  # Set individual probe_beta values to NA if the p-value is above max threshold
  #-----------------------------------------------------------------------------
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
    if (verbose) print(df_sum)
    verbosecat("\n")
  } else {
    verbosecat(sprintf("Probe filter: keeping all individual beta values that
                       exceed p-value.\n\n"))
  }
  # hist(rowSums(!is.na(filt_probe_beta_df))/ncol(filt_probe_beta_df),
  #      main = paste("Histogram of probe pass rate"),
  #      xlab = "Fraction passed probes")


  # Set entire row of probe measures to NA if not enough pass p-value check above
  #-----------------------------------------------------------------------------
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
  verbosecat(sprintf("Probe filter: discarding %.0f%% probes ( %i/ %i) because all measurements are now NA.
              %i probes now remain. \n\n",
      100*(nrow(filt_probe_beta_df) - nrow(filt_probe_beta_df2))/nrow(filt_probe_beta_df),
      nrow(filt_probe_beta_df) - nrow(filt_probe_beta_df2),
      nrow(filt_probe_beta_df),
      nrow(filt_probe_beta_df2)))


  # Discard probes with low design scores
  # df_design_score <- rbind(tdhia::design_scores$Pass_Score_Threshold,
  #                          tdhia::design_scores$Fail_Score_Threshold)
  # df_design_score <- tdhia::design_scores$Pass_Score_Threshold
  #
  # temp_filt_probe <- filt_probe_beta_df2 %>% tibble::rownames_to_column(var = "probe_id") %>%
  #   dplyr::mutate(probe_id = substr(.data$probe_id, 1, nchar(.data$probe_id)-1)) %>%
  # base::merge(y=dplyr::select(df_design_score, c("Assay_Design_Id", "Design_Score")),
  #       by.x = "probe_id", by.y = "Assay_Design_Id", all.x= TRUE, all.y= FALSE,
  #       no.dups = TRUE, sort = FALSE)
  #
  #
  # temp_filt_probe <- filt_probe_beta_df2 %>% tibble::rownames_to_column(var = "probe_id") %>%
  #   dplyr::left_join(y=dplyr::select(df_design_score, c("Assay_Design_Id", "Design_Score")),
  #         by =  dplyr::join_by(probe_id==Assay_Design_Id),
  #         keep = FALSE, na_matches = "never", #unmatched = "error",
  #         relationship = "one-to-one")


  # # Get design score for the remaining probes
  # a <-  filt_probe_beta_df2 %>% rownames_to_column(var = "probe_id") %>% select("probe_id")
  # a$probe_id <- substr(a$probe_id, 1, nchar(a$probe_id)-1)
  # b <- select(df_design_score, c("Assay_Design_Id", "Design_Score"))
  # c <- intersect(a$probe_id,b$Assay_Design_Id)
  #
  #
  # # get design score for all probes
  # a <-  tdhia::manifest_v1A2 %>% select("Probe_ID")
  # a$probe_id <- substr(a$Probe_ID, 1, nchar(a$Probe_ID)-1)
  # b <- select(df_design_score, c("Assay_Design_Id", "Design_Score"))
  # c <- intersect(a$probe_id,b$Assay_Design_Id)


  # keep_index <- temp_filt_probe$Design_Score >= min_probe_score
  # filt_probe_beta_df3 <- filt_probe_beta_df2[keep_index,]
  # filt_probe_pval_df3 <- filt_probe_pval_df2[keep_index,]
  # filt_design_scores3 <- temp_filt_probe$Design_Score[keep_index]
  # verbosecat(sprintf("Probe filter: discarding %.0f%% probes ( %i/ %i) b/c of
  #                    low design score < %.2f.\n\n",
  #                    100*sum(!keep_index)/nrow(filt_probe_beta_df3),
  #                    sum(!keep_index), nrow(filt_probe_beta_df3),
  #                    nrow(filt_probe_beta_df3)))

  # df_design_score$Assay_Design_Id


  # Report how many missing measurements at end of all filtering steps.
  verbosecat(sprintf("Probe filter: After all pruning and filtering, %.0f%% probes measurments ( %i/ %i) are now NA.\n\n",
                     100*sum(is.na(filt_probe_beta_df2)) / prod(dim(filt_probe_beta_df2)),
                     sum(is.na(filt_probe_beta_df2)),
                     prod(dim(filt_probe_beta_df2))))

 probe_beta <- list(probe_beta_df = filt_probe_beta_df2,
                    probe_pval_df = filt_probe_pval_df2,
                    # design_scores = filt_design_scores3,
                    platform = probe_beta$platform,
                    manifest = probe_beta$manifest,
                    probe_fail_rate_df = probe_fail_rate_df
  )

  return(probe_beta)


}



