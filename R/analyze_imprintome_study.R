




#' analyze_imprintome_study
#' @description performs a statistical analysis to identify CpG sites whose methylation
#' @param tbeta beta matrix
#' @param pheno study metadata
#'
#' @importFrom data.table :=
#'
#' @return sorted results from satistical analysis
#'
analyze_imprintome_study <- function (tbeta, pheno) {


  # Quantile Normalization



  # Beta matrix, rows:patients (samples), columns:cpg or ICR sites (features)
  tbetas<-tbeta
  dim(tbetas)
  dim(tbetas)
  # This line has to be true. This is the only thing you must care about when creating the two matrices.
  compare::compare(rownames(tbetas),rownames(pheno)) #TRUE


  # This is where you call
  start.time <- base::Sys.time()

  ind.hcc <-
    parallel::mclapply(
      stats::setNames(seq_len(ncol(tbetas)), dimnames(tbetas)[[2]]),
      GLMtest_2a,
      meth_matrix = tbetas,
      Y = pheno$case_control,
      X1 = pheno$sex,
      X2 = pheno$bw,
      mc.cores = 56
    )

  # End time
  end.time <- base::Sys.time()
  time.taken <- end.time - start.time

  # Setting metadata
  data.table::setattr(ind.hcc, 'class', 'data.frame')
  data.table::setattr(ind.hcc, "row.names", c(NA_integer_,4))
  data.table::setattr(ind.hcc, "names", make.names(names(ind.hcc), unique=TRUE))
  probelistnames <- names(ind.hcc)
  all.results <- t(data.table::as.data.table(ind.hcc))
  all.results <- data.table::as.data.table(all.results)
  all.results[, .data$probeID := probelistnames]
  stats::setNames(all.results, c("BETA","SE", "P_VAL", "probeID"))
  data.table::setcolorder(all.results, c("probeID","BETA","SE", "P_VAL"))

  # Organize and export resutls.
  all.results.sorted <- export_results(all.results, tbeta)

  return(all.results.sorted)

}



#' GLMtest_2a
#'
#' @description perform general linear modeling response variable Y to
#'
#' @param methcol numeric vector of which columns to use for predictors
#' @param meth_matrix matrix of predictor variables (beta values)
#' @param Y vector, response variable
#' @param X1 predictor ?
#' @param X2 predictor ?
#'
#'
#' @return fit of general linear model, including
#' - Estimate: point estimates of coefficients for each of the predictors
#' - Std. Error: standard error of the estimates
#' - Cumulative two-tailed probability
#'
GLMtest_2a = function(methcol, meth_matrix, Y, X1, X2 ) {

  mod = stats::glm(Y ~ meth_matrix[, methcol] + X1 + X2, family = "binomial")

  cf = summary(mod)$coefficients

  cf[2, c("Estimate", "Std. Error", "Pr(>|z|)")]

  return(cf)
}




#' export_results
#'
#' @description Organizes output from statistical analysis
#'
#' @param tbetas transposed beta_matrix
#' @param all.results output of general linear modeling analysis
#'
#' @importFrom rlang .data
#'
#' @return export
#'
#'
export_results <- function(all.results, tbetas) {

  #Add a column for the number of samples for each probe; this is just processing
  # the result and checking lambda to check if the assumption works ( generally,
  # anything close to 1 is good).
  # A Goodness-of-fit chi-square test is a statistical test used to determine
  # whether there is a significant difference between the observed and expected
  # frequencies in a categorical variable. It is commonly used in GLM (Generalized
  # Linear Model) analysis to assess the fit of the model to the data.

  # Transform methylation data again so that rows are probes and columns are samples
  tbetas_2<-t(tbetas)

  # Match order of all.results with order of probes in tbetas_1
  all.results<-all.results[match(rownames(tbetas_2),all.results$probeID),]
  all.results$N <- rowSums(!is.na(tbetas_2))
  # Convert results data table to Data.frame
  all.results <- base::as.data.frame(all.results) # High memory option
  # all.results <- Laurae::setDF(all.results)      # Low memory option, extra package req

  all.results.sorted<-all.results %>%
    dplyr::arrange(.data$P_VAL)

  utils::head(all.results.sorted)
  dplyr::count(all.results.sorted$P_VAL < 0.05)
  FileName<-paste("HCC_TruDx_v1_case_contorl_09152023.txt")

  #Lambda
  lambda <-
    stats::median(stats::qchisq( as.numeric( as.character( all.results.sorted$P_VAL)),
                        df = 1, lower.tail = F),
                        na.rm = T) / stats::qchisq( 0.5, 1)
  lambda # 1.56


  # export table of results as a .TXT.GZ FILE
  utils::write.table(all.results.sorted,file=FileName,na="NA", row.names = FALSE)

  return(all.results.sorted)
}








