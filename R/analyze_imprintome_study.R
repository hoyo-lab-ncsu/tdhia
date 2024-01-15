


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

  save(list = ls(all.names = TRUE), file = "analyze_imprintome_study.RData")
  # load(file = "analyze_imprintome_study.RData")


  # Quantile Normalization
  if (quantile_norm) {
    tbetas <- t(preprocessCore::normalize.quantiles(as.matrix(beta)))
  } else {
    tbetas <- t(as.matrix(beta))
  }
  rownames(tbetas) <- colnames(beta)
  colnames(tbetas) <- rownames(beta)

  # Reorder rows in tbetas so that ID matches order in rows of pheno (rownames for both)
  tbetas <- tbetas[order(match(rownames(tbetas),rownames(pheno))),]

  # Check that rownames between tbetas and pheno
  stopifnot(all(rownames(tbetas)==rownames(pheno)))




  # Threshold beta value to hemi-methylation
  hemi_tbetas <- (tbetas > 0.35) & (tbetas < 0.65)


  for_each_fun <-
    function (x) GLM_par(R = hemi_tbetas, Rind = x,   P = dplyr::select(pheno, cd_dry), Pind = 1,
                         C = dplyr::select(pheno, c(race_final, mat_bmi_lmp, smoking,
                                             smoke_preg)), family = "binomial")

  #Intiatialize parallel computing
  cl <- parallel::makePSOCKcluster(parallel::detectCores() - 1)
  doParallel::registerDoParallel(cl)

  start <- Sys.time()

  cat("Fitting model...\n")
  df_fits <- foreach::foreach(i=1:dim(hemi_tbetas)[2], .combine = rbind) %dopar% {
    for_each_fun(i)
  }
  colnames(df_fits) <-  c("Response", "Predictor", "BETA", "SE", "Z", "P_VAL",  "Formula")

  finish <- Sys.time()
  cat(sprintf("Processing time: "))
  finish - start


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


#' GLM_par
#'
#' @description perform general linear modeling with the equation of the form:
#'
#' Response ~ Predictor + Confounders
#' R[,Rind] ~ P[,Pind] + C[, 1] + C[, 2] + C[, 3] + C[, ...]
#'
#' One predictor and one response variable is assumed. Either R or P can be
#' parallelized (but only one of them)
#'
#' @param R dataframe where each column is a response variable.
#' @param Rind column index for R if it's a matrix (for parallel processing).
#' Then keep Rind = 1.
#' @param P dataframe where each column is a predictor variable.
#' @param Pind column index for P if it's a matrix (for parallel processing).
#' Then keep Pind = 1.
#' @param C matrix of confounder variables to be including in model
#' @param family string denoting GLM family
#'
#' @return fit of general linear model, including
#' - Estimate: point estimates of coefficients for each of the predictors
#' - Std. Error: standard error of the estimates
#' - Cumulative two-tailed probability
GLM_par = function(R, Rind = 1, P, Pind = 1, C, family = "binomial") {

  save(list = ls(all.names = TRUE), file = "GLM_par.RData")
  # load(file = "GLM_par.RData")

  stopifnot("GLM_par:error: R must be a dataframe" = is.data.frame(R) || is.matrix(R))
  stopifnot("GLM_par:error: P must be a dataframe" = is.data.frame(P) || is.matrix(P))
  if (length(Rind)!=1 || length(Pind)!=1) {
    stop(sprintf("GLM_par:error: both length() of Rind(==%.0f) and P(==%.0f) must be 1.",
                      length(Rind), length(Pind)))
  }

  # Formula string for confounders that preserve the column names (for record keeping)
  confounder_string <- paste0("C[, '",colnames(C),"']", c(rep(" +", ncol(C)-1), ""), collapse=" ")

  # Formula string for Response that preserve the column names (for record keeping)
  response_string <- paste0("R[, '", colnames(R)[Rind],"']")

  # Formula string for Predictor that preserve the column names (for record keeping)
  predictor_string <- paste0("P[, '", colnames(P)[Pind], "']")

  # Create a formula string dynamically using paste
  formula_string <- paste0(response_string, " ~ ", predictor_string, " + ", confounder_string)

  # Calculate glm
  mod = stats::glm(formula(formula_string), family = family)

  # Grab predictor from second row of summary output
  cf = summary(mod)$coefficients
  # Produce tidy dataframe of output
  out <- cbind(data.frame(Response = colnames(R)[Rind]),
               Predictor = colnames(P)[Pind],
               t(as.data.frame(cf[2,])),
               data.frame(Formula = formula_string))
  rownames(out) <- NULL


  return(out)
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








