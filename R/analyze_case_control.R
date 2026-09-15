


#' Compare Methylation Between Control and Case Groups
#'
#' Runs a two-sided Welch two-sample t-test for each methylation site and
#' summarizes group means, variability, missingness, and methylation changes.
#'
#' @param beta_matrix Numeric matrix or data frame with CpG or ICR sites in
#'   rows and samples in columns. Row names supply the site identifiers.
#' @param ctrl_cols Logical vector aligned to columns of beta_matrix.
#'   TRUE identifies controls (group 1); FALSE identifies cases (group 2).
#' @param n_adjust Number of comparisons passed to custom_p.adjust() with
#'   method = "fdr". Defaults to the original number of methylation sites.
#' @param filter_na_pval Logical; remove sites with missing raw p-values
#'   before multiple-testing adjustment. Defaults to TRUE.
#'
#' @details
#' A test is attempted only when each group has at least three nonmissing
#' measurements; otherwise its p-value is NA. Means and standard deviations
#' ignore missing values. Other t.test() errors, such as constant-data errors,
#' are not caught. Percent changes can be infinite or undefined when the
#' control mean is zero.
#'
#' @return A data frame in input-site order, optionally excluding missing
#'   p-values, with:
#'   - id and imp_id: original row number and site identifier.
#'   - p_val, adj_pval, adj_sig: raw p-value, FDR-adjusted p-value, and
#'     whether the adjusted p-value is below 0.05.
#'   - mean1, mean2, std1, std2: group means and standard deviations.
#'   - frac_na1 and frac_na2: missing fractions within each group.
#'   - delta: mean2 minus mean1.
#'   - perc_delta: 100 times delta divided by mean1.
#' @seealso [custom_p.adjust()]
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
