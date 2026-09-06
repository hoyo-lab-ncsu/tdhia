

#' convert_cpgs_to_icrs
#'
#' Converts a dataframe of beta values of specific CpG sites into beta
#' values for Imprint Control Region (ICR) sites. Calculates mean beta value for
#' cases where multiple CpG sites map to a single ICR site.
#'
#' @param cpg_beta a dataframe of beta values with rows representing cpg sites
#'  and columns representing different patients/ samples.
#' @param icr_mapping tba
#' @param sort_by_icr tba
#' @param max_icr_fail_rate if a fraction of samples for an ICR have NA values
#' above this proportion (0-1), the ICR is discarded from the results.
#' @param quantile_norm a boolean flag, when true normalizes the beta values
#' between columns of the output icr beta matrix. Default is FALSE because this
#' normalization is done earlier in the pipeline.
#' @param db_flag boolean when true export workspace to disk for debugging.
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @returns a dataframe containing beta values where rows are ICR ids and columns
#' are either (1) the basenames of the idat files or (2) some other mapping
#' specified when the data was loaded in \code{load_idat}.
#'
#' @export
convert_cpgs_to_icrs <- function(cpg_beta, icr_mapping = NULL, sort_by_icr = TRUE,
                                 max_icr_fail_rate = 0.20, quantile_norm = FALSE,
                                 db_flag = FALSE) {
  if (db_flag) {save(list = ls(all.names = TRUE), file = "convert_cpgs_to_icrs.RData")}
  # load(file = "convert_cpgs_to_icrs.RData")

  # Load dataframe that maps CpG sites to ICR site
  if (is.null(icr_mapping)) {icr_mapping = tdhia::mapping_cpg_icr_ids}

  # Look up the ICR id for each of the cpg sites, add it as a column
  #   This is used as a grouping variable for the next step.
  cpg_beta_df2 <- cpg_beta$cpg_beta_df %>%
      tibble::rownames_to_column( var = "CpG_id") %>%
      dplyr::left_join(,
                       y = dplyr::select(icr_mapping, c("ICR_id", "CpG_id")),
                       by = dplyr::join_by("CpG_id" == "CpG_id")) %>%
    dplyr::rename("ICR_ID" = "ICR_id")


  # Calculate mean beta value between cpg sites that map to same ICR site
  icr_beta_df <-
    cpg_beta_df2 %>%
    dplyr::select(-c("CpG_id")) %>%
    dplyr::group_by(.data$ICR_ID) %>%
    dplyr::summarize(dplyr::across(dplyr::where(is.numeric),
                                   function(x) mean(x, na.rm = TRUE)),
              n_CpGs = dplyr::n())


  # Sort ICR_ID by their number
  #   Extract icr id number and add as temp column, sort, remove temp column
  if (sort_by_icr) {
  icr_beta_df <- icr_beta_df %>%
    dplyr::mutate("icr_num_id" = as.numeric(gsub(".*_([0-9]+)$", "\\1", .data$ICR_ID))) %>%
    dplyr::arrange("icr_num_id") %>%
    dplyr::select(-c("icr_num_id"))
  }

  # Filter out entries that do not map to ICR region (icr_id == NA)
  icr_beta_df2 <- icr_beta_df[!is.na(icr_beta_df$ICR_ID),]

  # Filter ICRs that have too high of a fraction of failed measurements
  icr_discard <- rowSums(is.na(icr_beta_df2))/ncol(icr_beta_df2) > max_icr_fail_rate
  icr_beta_df3 <- icr_beta_df2[!icr_discard,]
  # Print out how many ICRs were discarded.
  cat(sprintf("ICR Filter: discarded %.0f%% of ICRs ( %i/ %i) b/c their signal fail rate was > %.f%%.
              %.0f ICRs still remain.\n",
              100*(nrow(icr_beta_df2)-nrow(icr_beta_df3))/length(unique(icr_mapping$ICR_id)),
              nrow(icr_beta_df2)-nrow(icr_beta_df3),
              length(unique(icr_mapping$ICR_id)),
              100*max_icr_fail_rate,
              nrow(icr_beta_df3)))

  icr_beta_df <- as.data.frame(icr_beta_df)

  # Remove entries that do not map to ICR
  icr_beta_df <- icr_beta_df[!is.na(icr_beta_df$ICR_ID),]

  # Convert ICR_ID to rownames
  rownames(icr_beta_df) <- icr_beta_df$ICR_ID
  # Extract number of CPGs for each ICR
  n_CpGs <- icr_beta_df$n_CpGs
  icr_beta_df <- icr_beta_df %>% dplyr::select(-c("ICR_ID", "n_CpGs"))

  icr_beta <- list(icr_beta_df = as.data.frame(icr_beta_df),
                   platform = cpg_beta$platform,
                   manifest = cpg_beta$manifest,
                   n_CpGs = n_CpGs,
                   cpg_icr_mapping = icr_mapping,
                   input_args = mget(setdiff(names(formals()), c("cpg_beta")), envir = environment()))

  return(icr_beta)
}



#
# # ICR Summary Table
# ################################################################################
# # Columns:
# #  icr_id: id number of ICR as annotated in the original epginetics paper
# #  mean_beta: mean of beta values across all patients
# #  std_beta: standard deviation of beta values across all patients
#
# summary_icr_beta_df <-
#   data.frame(icr_id = icr_beta_df$icr_id,
#              mean_beta = rowMeans(select(icr_beta_df, -icr_id, -n_CpGs)),
#              std_beta = rowSds(as.matrix(select(icr_beta_df, -icr_id, -n_CpGs),
#                                          na.rm=FALSE)))
# colnames(summary_icr_beta_df)<-
#   c(sprintf("ICR Name (n=%i)",ncol(select(icr_beta_df, -icr_id, -n_CpGs))),
#     "Mean (% Methylation)","Standard Deviation")
#
# write.csv(x = summary_icr_beta_df, file = paste0(output_dir_path, "/summary_icr_beta_df.csv"),
#           row.names = FALSE)
