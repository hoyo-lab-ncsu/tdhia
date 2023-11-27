

#' probe_beta_to_cpg
#'
#' Converts a dataframe of beta values of specific methylation probes into beta
#' values for CpG sites. Calculates mean beta value for cases where multiple
#' probes map to a single CpG site.
#'
#'
#' @param probe_beta_df a dataframe of beta values with rows representing cpg sites
#'  and columns representing different patients/ samples.
#' @param quantile_norm a boolean flag, when true normalizes the beta values
#' between columns of the icr beta dataframe. Default is FALSE because this
#' normalization is done earlier in the pipeline, at the probe level.
#'
#' @returns a dataframe containing beta values where rows are cpg sites and
#' columns are different patients/ samples.
#'
#' @export
#'
probe_beta_to_cpg <- function(probe_beta_df, quantile_norm = FALSE) {



  # Extract row names and add as first column
  probe_beta_matrix<-cbind(data.frame(probe_id = analyzed_probe_ids),probe_beta_matrix)
  rownames(probe_beta_matrix) <- NULL

  # Probe beta matrix
  ################################################################################
  # Write probe_beta_matrix to disk (all probes included, even nonexperiment ones)
  write.csv(x = probe_beta_matrix, file = paste0(output_dir_path, "/probe_beta_matrix.csv"),
            row.names = FALSE)




  # Preliminary setup for processing (1) CpG beta matrix and (2) ICR beta matrix
  ################################################################################
  # Load metadata file that maps CpG sites to ICR_id
  icr_df = read_csv(paste0(proj_path, "/reference_data/ICR_CpG_list.csv"))
  # Extract CpG_id from the CpG Probe labels, which is assumed to be all text
  # before the suffix
  #   Example: cg26679879_TC11" becomes "cg26679879"
  icr_df$CpG_id <- str_replace(icr_df$CpG_Probe,"_[a-zA-Z0-9]*","")

  # Add CpG site ID and ICR site ID columns to probe matrix for downstream processing
  temp_probe_beta_matrix <- probe_beta_matrix
  temp_probe_beta_matrix$cpg_id <- sapply(probe_beta_matrix$probe_id,
                                          function(x) mft$Name[which(x==mft$Probe_ID)[1]])




  # CpG Beta Matrix: average beta values between probe_id(s) that belong to same CpG site
  ################################################################################
  cpg_beta_matrix <-
    temp_probe_beta_matrix %>%
    # select(-icr_id) %>%
    group_by(cpg_id) %>%
    summarize(across(where(is.numeric),mean),
              n_probes=n())
  # Discard CpG IDs that do not map uniquely to a single genome location
  #   Define in manifest file where MAPINFO == 0, or NA
  f_unmmaped = function (x) is.na(x) | x==0;
  unmapped_cpg_ids <- sapply(cpg_beta_matrix$cpg_id,
                             function(x) f_unmmaped(mft$MAPINFO[which(x==mft$Name)[1]]))
  cpg_beta_matrix <- cpg_beta_matrix[!unmapped_cpg_ids,]
  # Write icr beta matrix to disk
  write.csv(x = cpg_beta_matrix, file = paste0(output_dir_path, "/cpg_beta_matrix.csv"),
            row.names = FALSE)





}
