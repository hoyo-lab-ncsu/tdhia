


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
  df <- data.frame(id = 1:nrow(beta_matrix), imp_id = NA, p_vals = rep(NA, nrow(beta_matrix)),
                   mean1 = NA, mean2 = NA, std1= NA, std2 = NA)
  

  #Loop through each row of beta_matrix
  for (i in 1:nrow(beta_matrix)) {
    #Split beta values into two groups based on ctrl_cols
    group1<- as.numeric(beta_matrix[i, ctrl_cols])
    group2<- as.numeric(beta_matrix[i, !ctrl_cols])

    #Perform two-sample t-test
    t_test_result<- stats::t.test(group1,group2)

    #Store the p-value in P_val
    df$p_vals[i] <- t_test_result$p.value
    df$mean1[i] <- mean(group1, na.rm = TRUE)
    df$mean2[i] <- mean(group2, na.rm = TRUE)
    df$std1[i] <-  sd(group1, na.rm = TRUE)
    df$std2[i] <-  sd(group2,  na.rm = TRUE)
    df$frac_na1[i] <-  sum(is.na(group1))/length(group1)
    df$frac_na2[i] <-  sum(is.na(group2))/length(group2)
  }


  # Adjust p-values based on number of ICRs (more liberal adjustment- this one has hits )
  df$adj_pval <- custom_p.adjust(df$p_vals, method = "fdr", n = n_adjust)
  df$adj_sig <- (df$adj_pval < 0.05)
  df$imp_id <- rownames(beta_matrix)
  df$delta <- df$mean2 - df$mean1
  df$perc_delta <- 100*(df$mean2 - df$mean1) / df$mean1
  # out <- data.frame(id = rownames(beta_matrix), p_val = p_vals,
  #                   adj_pval = adj_pvals, adj_sig = adj_sigs)
  return(df)
}
