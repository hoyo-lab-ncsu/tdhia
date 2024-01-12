


#' analyze_imprintome_study
#'
#' @description performs a statistical analysis with CpG sites whose methylation
#' state are correlated with a response variable while also considering several
#' covariates. CpG methylation status/ beta values are used as predictors, from
#' the study metadata a response variable is chosen along with any co-factors. The
#' presence of an association is tested for each of the predictor variables to
#' the response variable.
#'
#' @param beta a dataframe that is a beta matrix with CpG site IDs x patients
#' (row x col).
#' @param pheno a data frame of study metadata with patients x phenotype
#' variables (row x col). Contains the response variable Y and N covariates
#' X_1,...,X_N.
#' @param quantile_norm boolean flag, when true the beta matrix is quantile normalized
#' @importFrom data.table :=
#'
#' @return sorted results from statistical analysis
#'
analyze_imprintome_study <- function (beta, pheno, quantile_norm = TRUE) {


  # Quantile Normalization
  if (quantile_norm) {
    tbeta <- t(preprocessCore::normalize.quantiles(beta))
  } else {
    tbeta <- t(beta)
  }


  # Beta matrix, rows:patients (samples), columns:cpg or ICR sites (features)
  tbetas<-tbeta
  dim(tbetas)
  dim(tbetas)
  # This line has to be true. This is the only thing you must care about when creating the two matrices.
  stopifnot(compare::compare(rownames(tbetas),rownames(pheno)))


  # Perform general linear modeling of imprintome and study metadata
  start.time <- base::Sys.time()

  # Use all but 2 cores for parallel processing
  mc.cores <- max(c(parallel::detectCores()-2,1), na.rm = TRUE)

  ind.hcc <-
    parallel::mclapply(
      stats::setNames(seq_len(ncol(tbetas)), dimnames(tbetas)[[2]]),
      GLMtest_2a,
      meth_matrix = tbetas,
      Y = pheno$case_control,
      X1 = pheno$sex,
      X2 = pheno$bw,
      mc.cores = mc.cores
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
  all.results.sorted <- analyze_imprintome_export(all.results, tbeta)

  return(all.results.sorted)

}



#' GLMcall
#'
#' @description perform general linear modeling response variable Y to
#'
#' @param methcol vector of which column names to use for meth_matrix
#' @param meth_matrix matrix of predictor variables (beta values)
#' @param Y vector, response variable
#' @param ... any additional covariates variables
#'
#'
#' @return fit of general linear model, including
#' - Estimate: point estimates of coefficients for each of the predictors
#' - Std. Error: standard error of the estimates
#' - Cumulative two-tailed probability
#'
GLMcall = function(methcol, meth_matrix, Y, ... ) {
  args = list(...)
  # Create a formula string dynamically using paste
  formula_string <- paste("Y ~ meth_matrix[, methcol]", paste(args, collapse = " + "), sep = " + ")

  mod = stats::glm(formula_string, family = "binomial")

  cf = summary(mod)$coefficients

  cf[2, c("Estimate", "Std. Error", "Pr(>|z|)")]

  return(cf)
}




#' analyze_imprintome_export
#'
#' @description Organizes output from statistical analysis for export
#'
#' @param tbetas transposed beta_matrix
#' @param all.results output of general linear modeling analysis
#'
#' @importFrom rlang .data
#'
#' @return export
#'
#'
analyze_imprintome_export <- function(all.results, tbetas) {

  #Add a column for the number of samples for each probe; this is just processing
  # the result and checking lambda to check if the assumption works ( generally,
  # anything close to 1 is good).
  # A Goodness-of-fit chi-square test is a statistical test used to determine
  # whether there is a significant difference between the observed and expected
  # frequencies in a categorical variable. It is commonly used in GLM (Generalized
  # Linear Model) analysis to assess the fit of the model to the data.

  # Transform methylation data again so to probes x samples (row x col)
  tbetas_2<-t(tbetas)

  # Match order of all.results with order of probes in tbetas_1
  all.results<-all.results[match(rownames(tbetas_2),all.results$probeID),]
  all.results$N <- rowSums(!is.na(tbetas_2))
  # Convert results data table to data.frame
  all.results <- base::as.data.frame(all.results) # High memory option
  # all.results <- Laurae::setDF(all.results)      # Low memory option, extra package req

  # Sort rows by p value
  all.results.sorted<-all.results %>%
    dplyr::arrange(.data$P_VAL)

  utils::head(all.results.sorted)
  dplyr::count(all.results.sorted$P_VAL < 0.05)
  FileName<-paste("imprintome_study_results.txt")

  # Lambda
  lambda <-
    stats::median(stats::qchisq( as.numeric( as.character( all.results.sorted$P_VAL)),
                        df = 1, lower.tail = F),
                        na.rm = T) / stats::qchisq( 0.5, 1)
  lambda # 1.56


  # export table of results as a .TXT.GZ FILE
  utils::write.table(all.results.sorted,file=FileName,na="NA", row.names = FALSE)

  return(all.results.sorted)
}








