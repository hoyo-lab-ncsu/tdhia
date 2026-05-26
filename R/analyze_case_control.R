


#' analyze_case_control
#'
#' @description analyze for changes in beta value between a control group and an
#' experiment group.
#'
#' @param beta_matrix a matrix of beta values, cpg/ ICR site (rows) by patients (columns)
#' @param ctrl_cols boolean vector of which columns in beta matrix are control group
#' @param n_adjust number to adjust p-values for FDR multiple comparisons
#' @param filter_na_pval todo
#' @export
analyze_case_control <- function(beta_matrix, ctrl_cols, n_adjust = nrow(beta_matrix), filter_na_pval = TRUE) {
  # Initialize a vector of NAs to store p-values
  df <- data.frame(id = 1:nrow(beta_matrix), imp_id = NA, p_val = rep(NA, nrow(beta_matrix)),
                   adj_pval = NA, adj_sig = NA, mean1 = NA, mean2 = NA, std1= NA, std2 = NA)
  

  #Loop through each row of beta_matrix
  for (i in 1:nrow(beta_matrix)) {
    #Split beta values into two groups based on ctrl_cols
    group1<- as.numeric(beta_matrix[i, ctrl_cols])
    group2<- as.numeric(beta_matrix[i, !ctrl_cols])

    #Perform two-sample t-test
    if ((length(group1[!is.na(group1)]) > 2) & 
        (length(group2[!is.na(group2)]) > 2)   ) {
      t_test_result<- stats::t.test(group1,group2)
    } else {
      t_test_result = list(p.value = NA)
    }

    #Store the p-value in P_val
    df$p_val[i] <- t_test_result$p.value
    df$mean1[i] <- mean(group1, na.rm = TRUE)
    df$mean2[i] <- mean(group2, na.rm = TRUE)
    df$std1[i] <-  stats::sd(group1, na.rm = TRUE)
    df$std2[i] <-  stats::sd(group2,  na.rm = TRUE)
    df$frac_na1[i] <-  sum(is.na(group1))/length(group1)
    df$frac_na2[i] <-  sum(is.na(group2))/length(group2)
  }

  df$delta <- df$mean2 - df$mean1
  df$perc_delta <- 100*(df$mean2 - df$mean1) / df$mean1
  df$imp_id <- rownames(beta_matrix)
  
  # Remove rows where a p-value can't be calculated (required for p-value adjustment)
  if (filter_na_pval) {    df <- df[!is.na(df$p_val),]    }
  # Adjust p-values based on number of ICRs (more liberal adjustment- this one has hits )
  df$adj_pval <- custom_p.adjust(df$p_val, method = "fdr", n = n_adjust)
  df$adj_sig <- (df$adj_pval < 0.05)
  

  
  return(df)
}
