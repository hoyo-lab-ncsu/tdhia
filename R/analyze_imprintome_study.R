


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
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
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

  # Helper function for parallel processing of fitting models
  for_each_fun <-
    function (x) GLM_parallel(R = hemi_tbetas, Rind = x,   P = dplyr::select(pheno, cd_dry), Pind = 1,
                         C = dplyr::select(pheno, c(race_final, mat_bmi_lmp, smoking,
                                             smoke_preg)), family = "binomial")

  #Initialize parallel computing
  cl <- parallel::makePSOCKcluster(parallel::detectCores() - 1)
  doParallel::registerDoParallel(cl)

  # Fit mdoels to data in parallel
  start <- Sys.time()
  cat("Fitting model...\n")
  df_fits <- foreach::foreach(i=1:dim(hemi_tbetas)[2], .combine = rbind) %dopar% {
    for_each_fun(i)
  }
  colnames(df_fits) <-  c("Response", "Predictor", "EST", "SE", "Z", "P_VAL",  "Formula")
  finish <- Sys.time()
  cat(sprintf("Processing time: %f %s\n", finish - start, units(finish - start)))

  # Sort by p value
  df_fits_sorted <- df_fits %>% dplyr::arrange(.data$P_VAL)

  # Adjust pvalues for FDR correction
  df_fits_sorted$ADJ_P_VAL <- stats::p.adjust(df_fits_sorted$P_VAL, method = "fdr")

  cat(sprintf("%.0f cpg sites have p_val < 0.05\n", sum(df_fits_sorted$P_VAL < 0.05)))
  cat(sprintf("%.0f cpg sites have adj_p_val < 0.05\n", sum(df_fits_sorted$ADJ_P_VAL < 0.05)))

  # TODO: fill in what this does
  lambda <-
    stats::median(stats::qchisq( as.numeric( as.character(df_fits_sorted$P_VAL)),
                                 df = 1, lower.tail = F),
                  na.rm = T) / stats::qchisq( 0.5, 1)
  cat(sprintf("Lambda: %f\n", lambda))


  return(df_fits_sorted)
}



#' GLM_parallel
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
GLM_parallel = function(R, Rind = 1, P, Pind = 1, C, family = "binomial") {

  save(list = ls(all.names = TRUE), file = "GLM_parallel.RData")
  # load(file = "GLM_parallel.RData")

  stopifnot("GLM_parallel:error: R must be a dataframe" = is.data.frame(R) || is.matrix(R))
  stopifnot("GLM_parallel:error: P must be a dataframe" = is.data.frame(P) || is.matrix(P))
  if (length(Rind)!=1 || length(Pind)!=1) {
    stop(sprintf("GLM_parallel:error: both length() of Rind(==%.0f) and P(==%.0f) must be 1.",
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




