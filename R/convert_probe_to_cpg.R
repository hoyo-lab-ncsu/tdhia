

#' convert_probe_to_cpg
#'
#' Converts a dataframe of beta values of specific methylation probes into beta
#' values for CpG sites. Calculates mean beta value for cases where multiple
#' probes map to a single CpG site.
#'
#' @param probe_beta_df a dataframe of beta values with rows representing cpg sites
#'  and columns representing different patients/ samples.
#' @param quantile_norm a boolean flag, when true normalizes the beta values
#' between columns of the icr beta dataframe. Default is FALSE because this
#' normalization is done earlier in the pipeline, at the probe level.
#' @param mft manifest file for particular imprintome array used in data collection.
#' A dataframe with metadata mapping probe_ids to CpG sites. See ?manifest_v1A2 dataset.
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
convert_probe_to_cpg <- function(probe_beta_df, mft = NULL, quantile_norm = FALSE,
                                 discard_unmapped_cpgs = TRUE) {

  # # Convert rownames to prode_beta_df to the first column
  # probe_beta_df2<-cbind(data.frame(Probe_ID = rownames(probe_beta_df)),probe_beta_df)
  # rownames(probe_beta_df2) <- NULL

  if (is.null(mft)) {mft = tdhia::manifest_v1A2}

  # Match CpG site ID for each probe ID in probe_beta_df2
  probe_beta_df2 <-
    dplyr::select(
      dplyr::left_join(probe_beta_df,
                       dplyr::select(mft, c("Probe_ID", "Name")),
                       "Probe_ID",-"Probe_ID"), -"Probe_ID"
      ) %>%
    dplyr::rename("CpG_ID" = "Name")

  # CpG Beta Matrix: average beta values between probe_id(s) that belong same CpG site
  cpg_beta_df <-
    probe_beta_df2 %>%
    dplyr::group_by(.data$CpG_ID) %>%
    dplyr::summarize(dplyr::across(dplyr::where(is.numeric),mean),
              n_probes = dplyr::n())

  # Discard CpG_ID(s) that do not map uniquely to a single genome location
  #   Defined in manifest file where MAPINFO == (0 or NA)
  if (discard_unmapped_cpgs) {
  f_unmmaped = function (x) is.na(x) | x==0;
  unmapped_cpg_ids <- sapply(cpg_beta_df$CpG_ID,
                             function(x) f_unmmaped(mft$MAPINFO[which(x==mft$Name)[1]]))
  cpg_beta_df <- cpg_beta_df[!unmapped_cpg_ids,]
  }

  return(cpg_beta_df)
}
