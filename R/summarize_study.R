
#' summarize_study
#'
#' @description Prints out results of a study analysis from analyze_association.
#' @param dfs a list of dataframes that give model output results for each
#' variable in model.
#' @param varnames vector of strings specifying which variables to cummarize.
#' @param max_p_val maximum p-value threshold for reported results. Both
#' adjusted and unadjusted p-values are reported.
#' @param print_sites boolean, when true, will print all cpg/icr sites that are
#' statistically significant (adjust p-value)
#' @param print_confounders boolean, when true, will print out statistical
#' summary of confoudner variables
#' @export
summarize_study <- function(dfs, varnames = NULL, max_p_val = 0.05,
                            print_sites = TRUE, print_confounders = FALSE) {
  cat(sprintf("Formula: %s \n", dfs$example_formula))

  if (is.null(varnames)) varnames <- names(dfs)

  sig_list = list()
  for (n in seq_along(varnames)) {
    if (dfs[[n]]$Confounder[1] == 0 || print_confounders ) {
      cat(sprintf("%s:\n", varnames[n]))
      cat(sprintf(">>  %.0f imprint sites have p_val < %.2f\n",
                  sum(dfs[[n]]$P_VAL < max_p_val), max_p_val))
      cat(sprintf(">>  %.0f imprint sites have adj_p_val < %.2f\n",
                  sum(dfs[[n]]$ADJ_P_VAL < max_p_val), max_p_val))
      if (sum(dfs[[n]]$ADJ_P_VAL < max_p_val)>0 && print_sites) {

        print(dfs[[n]][dfs[[n]]$ADJ_P_VAL < max_p_val,] %>%
                dplyr::select(,-c("Formula", "Model_Id")))
      }

      # Print ho wmany model fittings failed
      cat(sprintf("%.0f/ %.0f of model fits failed.\n", sum(is.na(dfs[[n]]$Estimate)), 
                  nrow(dfs[[n]])))
      cat("\n")
    }

    sig_list[[n]] <- dfs[[n]][ dfs[[n]]$ADJ_P_VAL < max_p_val, ]
  }

  df_sig = do.call(rbind, sig_list)
}

