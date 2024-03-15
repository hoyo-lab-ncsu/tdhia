
#' summarize_study
#'
#' @description
#'  Prints out results of a study analysis.
#' @param dfs a list of dataframes that give model output results for each
#' variable in model.
#' @param max_p_val maximum p-value threshold for reported results. Both
#' adjusted and unadjusted p-values are reported.
#' @param print_sites boolean, when true, will print all cpg/icr sites that are
#' statistically significant (adjust p-value)
#' @param print_confounders boolean, when true, will print out statistical
#' summary of confoudner variables
summarize_study <- function(dfs, max_p_val, print_sites = TRUE, print_confounders = FALSE) {
  cat("Summarizing Results...\n")
  model_vars <- names(dfs)
  for (n in seq_along(model_vars)) {
    if (dfs[[n]]$Confounder[1] == 0 || print_confounders) {
      cat(sprintf("%s:\n", model_vars[n]))
      cat(sprintf(">>  %.0f cpg sites have p_val < %.2f\n",
                  sum(dfs[[n]]$P_VAL < max_p_val), max_p_val))
      cat(sprintf(">>  %.0f cpg sites have adj_p_val < %.2f\n",
                  sum(dfs[[n]]$ADJ_P_VAL < max_p_val), max_p_val))
      if (sum(dfs[[n]]$ADJ_P_VAL < max_p_val)>0 && print_sites) {

        print(dfs[[n]][dfs[[n]]$ADJ_P_VAL < max_p_val,] %>%
                dplyr::select(,-c("Formula", "Model_Id")))
      }
      cat("\n")
    }
  }
}

