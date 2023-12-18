

#' convert_probes_to_cpgs
#'
#' Converts a dataframe of beta values of specific methylation probes into beta
#' values for CpG sites. Calculates mean beta value for cases where multiple
#' probes map to a single CpG site.
#'
#' @param probe_beta a list object that must contain
#'        - probe_beta_df: beta value dataframe, probe_id x patient_id
#'        - probe_pval_df: p-value value dataframe, probe_id x patient_id
#'        - Probe_ID: vector of string containing ID for each probe
#' @param quantile_norm a boolean flag, when true normalizes the beta values
#' between columns of the icr beta dataframe. Default is FALSE because this
#' normalization is done earlier in the pipeline, at the probe level.
#' @param discard_unmapped_cpgs  a boolean flag, when set to TRUE, discards any
#' CpG sites that are not mapped to a unique genomic location from MAPINFO column
#' in the manifest file.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @returns a dataframe containing beta values where rows are cpg sites and columns
#' are either (1) the basenames of the idat files or (2) some other mapping
#' specified when the data was loaded in \code{load_idat}.
#'
#'
convert_probes_to_cpgs <- function(probe_beta, quantile_norm = FALSE,
                                 discard_unmapped_cpgs = TRUE) {
  # Get manifest data
  mft = probe_beta$manifest


  # Match CpG site ID for each probe ID in probe_beta_df2
  probe_beta_df2 <-
    probe_beta$probe_beta_df %>%
    tibble::rownames_to_column(var = "Probe_ID") %>%
    dplyr::left_join(y = dplyr::select(mft, c("Probe_ID", "Name")),
                     "Probe_ID",-"Probe_ID") %>%
    dplyr::select(-c("Probe_ID")) %>%
    dplyr::rename("CpG_ID" = "Name")


  # CpG Beta Matrix: average beta values between probe_id(s) that belong same CpG site
  cpg_beta_df <-
    probe_beta_df2 %>%
    dplyr::group_by(.data$CpG_ID) %>%
    dplyr::summarize(dplyr::across(dplyr::where(is.numeric),
                                   function(x) mean(x, na.rm = TRUE)),
                     n_probes = dplyr::n())
  # Report how many CpGs discarded from manifest
  n_unique_cpgs_mft <- length(unique(mft$Name))
  n_unique_cpgs_data <- length(unique(cpg_beta_df$CpG_ID))
  cat(sprintf("CpG filter: discarded %.0f%% of CpG entries ( %i/ %i) from manifest because
              they have no associated probe data that passed probe filtering.
              %i CpG entries remain.\n",
              100*(n_unique_cpgs_mft-n_unique_cpgs_data)/n_unique_cpgs_mft,
              n_unique_cpgs_mft-n_unique_cpgs_data, n_unique_cpgs_data,
              n_unique_cpgs_data))



  # Discard CpG_ID(s) that do not map uniquely to a single genome location
  #   Defined in manifest file where MAPINFO == (0 or NA)
  if (discard_unmapped_cpgs) {
  f_unmmaped = function (x) is.na(x) | x==0;
  unmapped_cpg_ids <- sapply(cpg_beta_df$CpG_ID,
                             function(x) f_unmmaped(mft$MAPINFO[which(x==mft$Name)[1]]))

  cpg_beta_df <- cpg_beta_df[!unmapped_cpg_ids,]
  cat(sprintf("CpG filter: discarded %.0f%% of CpG entries (%i/ %i) because they
      do not map uniquely to the genome. %i CpG sites remain.\n",
      100*sum(unmapped_cpg_ids)/length(unmapped_cpg_ids), sum(unmapped_cpg_ids),
      length(unmapped_cpg_ids), length(unmapped_cpg_ids) - sum(unmapped_cpg_ids)))

  }

  cpg_beta <- list(cpg_beta_df = cpg_beta_df,
                     platform = probe_beta$platform,
                     manifest = probe_beta$manifest)

  return(cpg_beta)
}
