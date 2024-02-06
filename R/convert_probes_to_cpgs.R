

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
#' @param discard_non_icr_cpgs a boolean flag, when set to TRUE, discards any CpG
#'  sites that do not map to an ICR.
#' @param return_n_probes a boolean flag, when set to TRUE, adds a column,n_probes,
#'  to the output cpg beta matrix (cpg_beta_df), that represents the
#'  number of probes whose beta value was averaged together for that cpg site.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @returns a dataframe containing beta values where rows are cpg sites and columns
#' are either (1) the basenames of the idat files or (2) some other mapping
#' specified when the data was loaded in \code{load_idat}.
#'
convert_probes_to_cpgs <- function(probe_beta, quantile_norm = FALSE,
                                   discard_unmapped_cpgs = TRUE,
                                   discard_non_icr_cpgs = TRUE,
                                   return_n_probes = FALSE) {
  # Get manifest data
  mft = probe_beta$manifest
  icr_mapping = tdhia::mapping_cpg_icr_ids

  save(list = ls(all.names = TRUE), file = "convert_probes_to_cpgs.RData")
  # load(file = "convert_probes_to_cpgs.RData")


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
  cat(sprintf("CpG filter: %.0f of probes in dataset mapped to %.0f unique CpG sites id data.\n",
              nrow(probe_beta$probe_beta_df), nrow(cpg_beta_df)))

  # Convert to dataframe to add rownames and drop first column (Cpg IDs)
  cpg_beta_df <- as.data.frame(cpg_beta_df)
  rownames(cpg_beta_df) <- cpg_beta_df$CpG_ID
  cpg_beta_df <- dplyr::select(cpg_beta_df,-.data$CpG_ID)


  # Discard CpG_ID(s) that do not map uniquely to a single genome location
  #   Defined in manifest file where MAPINFO == (0 or NA)
  if (discard_unmapped_cpgs) {
    f_unmmaped = function (x) is.na(x) | x==0;
    unmapped_cpg_ids <- sapply(rownames(cpg_beta_df),
                               function(x) f_unmmaped(mft$MAPINFO[which(x==mft$Name)[1]]))

    # Remove unmapped cpgs
    cpg_beta_df <- cpg_beta_df[!unmapped_cpg_ids,]

    cat(sprintf("CpG filter: discarded %.0f%% of CpG sites (%i/ %i) because they
      do not map uniquely to the genome. %i CpG sites remain.\n",
                100*sum(unmapped_cpg_ids)/length(unmapped_cpg_ids), sum(unmapped_cpg_ids),
                length(unmapped_cpg_ids), nrow(cpg_beta_df)))

  }


  if (discard_non_icr_cpgs){
    discard_non_icr <- cpg_beta_df %>%
      tibble::rownames_to_column("CpG_id") %>%
      dplyr::left_join(y = dplyr::select(icr_mapping, c("CpG_id", "ICR_id")),
                              by = "CpG_id",unmatched = "drop",
                multiple = "first") %>%
      dplyr::pull("ICR_id") %>% is.na()
    # Remove non-icr cpgs
    cpg_beta_df <- cpg_beta_df[!discard_non_icr,]

    cat(sprintf("CpG filter: discarded %.0f%% of CpG sites (%i/ %i) because they
      do not map to an ICR. %i CpG sites remain.\n",
                100*sum(discard_non_icr)/nrow(cpg_beta_df), sum(discard_non_icr),
                nrow(cpg_beta_df), nrow(cpg_beta_df)))
  }


  # Extract n_probes from dataframe (to be stored separately)
  n_probes <- cpg_beta_df$n_probes
  cpg_beta_df <- cpg_beta_df %>% dplyr::select(-n_probes)

  if (quantile_norm) {
    # Quantile Normalization by column
    cat("Performing quantile normalization...\n")
    norm_mat <- preprocessCore::normalize.quantiles(as.matrix(cpg_beta_df))
    dimnames(norm_mat) <- list(rownames(cpg_beta_df), colnames(cpg_beta_df))
    cpg_beta_df <- norm_mat
  }


  cpg_beta <- list(cpg_beta_df = as.data.frame(cpg_beta_df),
                   platform = probe_beta$platform,
                   manifest = probe_beta$manifest,
                   cpg_n_probes = n_probes)

  return(cpg_beta)
}
