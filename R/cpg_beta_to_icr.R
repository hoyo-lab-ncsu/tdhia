

#' cpg_beta_to_icr
#'
#' Converts a dataframe of beta values of specific CpG sites into beta
#' values for Imprint Control Region (ICR) sites. Calculates mean beta value for
#' cases where multiple CpG sites map to a single ICR site.
#'
#' @param cpg_beta_df a dataframe of beta values with rows representing cpg sites
#'  and columns representing different patients/ samples.
#' @param quantile_norm a boolean flag, when true normalizes the beta values
#' between columns of the icr beta dataframe. Default is FALSE because this
#' normalization is done earlier in the pipeline, at the probe level.
#'
#' @returns a dataframe containing beta values where rows are probe_ids and columns
#' are the basenames of the idat files or some other mapping specifed by column_mappings.
#'
#' @export
#'
cpg_beta_to_icr <- function(cpg_beta_df, quantile_norm = FALSE) {



  # ICR Beta Matrix: average beta values between cpg_id(s) that belong to same ICR site
  ################################################################################
  # Warning: do not average all probe betas that belong to an ICR site
  #   Instead, average all CPG betas that belong to an ICR site

  # Look up the ICR id for each of the cpg sites, add it as a column
  #   This is used as a grouping variable for the next step.
  temp_cpg_beta_matrix <- cpg_beta_matrix
  temp_cpg_beta_matrix$icr_id <- unname(sapply(cpg_beta_matrix$cpg_id,
                                               function(x) icr_df$ICR_id[which(x==icr_df$CpG_id)[1]]))
  icr_beta_matrix <-
    temp_cpg_beta_matrix %>%
    select(-cpg_id, -n_probes) %>%
    group_by(icr_id) %>%
    summarize(across(where(is.numeric),mean),
              n_CpGs=n())

  # Sort ICR_ID by their number
  #   Extract icr id number and add as temp column, sort, remove temp column
  icr_beta_matrix <- icr_beta_matrix %>%
    mutate(icr_num_id = as.numeric(gsub(".*_([0-9]+)$", "\\1", icr_id))) %>%
    arrange(icr_num_id) %>%
    select(-icr_num_id)

  # Filter out entries that do not map to ICR region (icr_id == NA)
  icr_beta_matrix <- icr_beta_matrix[!is.na(icr_beta_matrix$icr_id),]
  # Write icr beta matrix to disk
  write.csv(x = icr_beta_matrix, file = paste0(output_dir_path, "/icr_beta_matrix.csv"),
            row.names = FALSE)





}



#
# # ICR Summary Table
# ################################################################################
# # Columns:
# #  icr_id: id number of ICR as annotated in the original epginetics paper
# #  mean_beta: mean of beta values across all patients
# #  std_beta: standard deviation of beta values across all patients
#
# summary_icr_beta_matrix <-
#   data.frame(icr_id = icr_beta_matrix$icr_id,
#              mean_beta = rowMeans(select(icr_beta_matrix, -icr_id, -n_CpGs)),
#              std_beta = rowSds(as.matrix(select(icr_beta_matrix, -icr_id, -n_CpGs),
#                                          na.rm=FALSE)))
# colnames(summary_icr_beta_matrix)<-
#   c(sprintf("ICR Name (n=%i)",ncol(select(icr_beta_matrix, -icr_id, -n_CpGs))),
#     "Mean (% Methylation)","Standard Deviation")
#
# write.csv(x = summary_icr_beta_matrix, file = paste0(output_dir_path, "/summary_icr_beta_matrix.csv"),
#           row.names = FALSE)
