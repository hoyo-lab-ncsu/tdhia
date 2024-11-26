


#' analyze_case_control
#'
#' @description analyze for changes in beta value between a control group and an
#' experiment group.
#'
#' @param beta_matrix a matrix of beta values, cpg/ ICR site (rows) by patients (columns)
#' @param ctrl_cols boolean vector of which columns in beta matrix are control group
#' @param n_adjust number to adjust p-values for FDR multiple comparisons
#'
#' @export
analyze_case_control <- function(beta_matrix, ctrl_cols, n_adjust = nrow(beta_matrix)) {
  #Initialize a vector of NAs to store p-values
  p_vals <- rep(NA, nrow(beta_matrix))

  #Loop through each row of beta_matrix
  for (i in 1:nrow(beta_matrix)) {
    #Split beta values into two groups based on ctrl_cols
    group1<- beta_matrix[i, ctrl_cols]
    group2<- beta_matrix[i, !ctrl_cols]

    #Perform two-sample t-test
    t_test_result<- stats::t.test(group1,group2)

    #Store the p-value in P_val
    p_vals[i] <- t_test_result$p.value
  }


  # Adjust p-values based on number of ICRs (more liberal adjustment- this one has hits )
  adj_pvals <- custom_p.adjust(p_vals, method = "fdr", n = n_adjust)
  adj_sigs <- (adj_pvals < 0.05)
  out <- data.frame(id = rownames(beta_matrix), p_val = p_vals,
                    adj_pval = adj_pvals, adj_sig = adj_sigs)
  return(out)
}
