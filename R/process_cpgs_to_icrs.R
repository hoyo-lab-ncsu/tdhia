

#' Aggregate CpG Methylation to ICRs
#'
#' Maps CpG sites to imprinting control regions and averages their beta
#' values separately for each sample, ignoring missing values.
#'
#' @param cpg_beta List returned by convert_probes_to_cpgs(), containing
#'   cpg_beta_df (CpGs in rows, samples in columns), platform, and manifest.
#' @param icr_mapping Data frame containing CpG_id and ICR_id columns.
#'   NULL uses the package's mapping_cpg_icr_ids data.
#' @param sort_by_icr Logical; request sorting by the numeric ICR suffix.
#'   The current implementation passes a quoted column name to arrange(),
#'   so the requested numeric ordering is not reliably applied.
#' @param max_icr_fail_rate Numeric missing-value fraction used to calculate
#'   and print an ICR discard summary. The filtered intermediate table is
#'   currently not used to construct the return value, so this argument does
#'   not remove ICRs from the returned data.
#' @param quantile_norm Logical; currently unused. No quantile normalization
#'   is performed, regardless of this argument.
#' @param db_flag Logical; save the initial environment to
#'   convert_cpgs_to_icrs.RData in the working directory.
#'
#' @details
#' CpGs are joined to the supplied mapping before averaging; multiple mapping
#' rows can therefore contribute multiple times. Unmapped ICR entries are
#' removed. An all-missing set of beta values produces NaN. The reported
#' failure fraction is calculated before dropping the ICR ID and count columns.
#' A filtering summary is printed on every call.
#'
#' @return A named list with:
#'   - icr_beta_df: data frame of mean beta values, ICR IDs in rows and
#'     sample IDs in columns.
#'   - platform and manifest: metadata copied from cpg_beta.
#'   - n_CpGs: number of joined CpG rows per returned ICR, in table order;
#'     this is not a per-sample count of nonmissing measurements.
#'   - cpg_icr_mapping: mapping table used for aggregation.
#'   - input_args: argument values excluding cpg_beta.
#' @seealso [convert_probes_to_cpgs()], [tdhia_pipeline()]
#' @importFrom magrittr %>%
#' @importFrom rlang .data
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
