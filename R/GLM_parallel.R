
#' GLM_parallel
#'
#' @description perform general linear modeling with the equation of the form:
#'
#' Response ~ Predictor + Confounders
#' \code{R[,Rind] ~ P[,Pind] + C[, 1] + C[, 2] + C[, 3] + C[, ...]}
#'
#' One predictor and one response variable (column) is assumed. If fitting a
#' series of models for several R's or P's, this function can be called in
#' parallel, and the particular column to be used for each iteration can be
#' specified by the column index Rind or Pind. NA values from R or P can be
#' imputed with the impute_na flag.
#'
#' @param R dataframe where each column is a Response variable.
#' @param Rind column index for R if it has multiple columns (for parallel processing).
#' Default for Rind is 1 for single column dataframe.
#' @param P dataframe where each column is a Predictor variable.
#' @param Pind column index for P if it has multiple columns (for parallel processing).
#' Default for Pind is 1 for single column dataframe.
#' @param C dataframe of confounder variable(s) to be including in model
#' @param family string denoting GLM family
#' @param verbose boolean flag, when true prints fits to model.
#' @param impute_na boolean flag, when TRUE inputs missing NA values with MICE
#' package.
#' @param db_flag boolean flag, when TRUE saves workspace to disk for debugging
#'
#' @return fit of general linear model, including
#' - Estimate: point estimates of coefficients for each of the predictors
#' - Std. Error: standard error of the estimates
#' - Cumulative two-tailed probability
GLM_parallel = function(R, Rind = 1, P, Pind = 1, C, family = "binomial",
                        verbose = FALSE, impute_na = TRUE, db_flag = FALSE) {
  if (db_flag) {save(list = ls(all.names = TRUE), file = "GLM_parallel.RData") }
  # load(file = "GLM_parallel.RData")

  stopifnot("GLM_parallel:error: R must be a dataframe" = is.data.frame(R) ||
              is.matrix(R))
  stopifnot("GLM_parallel:error: P must be a dataframe" = is.data.frame(P) ||
              is.matrix(P))
  if (length(Rind)!=1 || length(Pind)!=1) {
    stop(sprintf(paste0("GLM_parallel:error: both length() of Rind(==%.0f) and ",
                        "Pind(==%.0f) must be 1."),
                 length(Rind), length(Pind)))
  }

  # Condense R and P and C into a single matrix: R,P,C
  model_data <- cbind(data.frame(R = R[, Rind], P = P[, Pind]), C)
  # Rename Response and Predictor columns to their variables
  names(model_data)[1] <- colnames(R)[Rind]
  names(model_data)[2] <- colnames(P)[Pind]

  # model_data[runif(10,1,nrow(model_data)),1]<- NA


  # Calculate the number of imputations required
  # Generally, number of imputations is % of missing data
  # Unless no missing data, than imputes is set to 1 to not error package
  # If n_imputes <5, then set to 5, because that is considered the min if there
  # are any missing values
  n_imputes <- ceiling(100 * sum(is.na( model_data[,1:2])) /
                         (2 * nrow(model_data)))
  # If imputation is needed, set to at least 5 (original paper rec)
  if(n_imputes > 0 && n_imputes < 5) {n_imputes = 5}

  # Define formula string for fitted model
  formula_string <-
    paste0(colnames(R)[Rind], " ~ ", colnames(P)[Pind]," + ",
           paste0(colnames(C), c(rep(" +", ncol(C)-1), ""), collapse=" "))

  # Impute missing data if flag is set and data actually missing
  if (impute_na && n_imputes > 0 ) {

    # Perform multiple imputations of the dataset
    imp <- mice::mice(model_data, print = FALSE, m = n_imputes, maxit = 10, seed = 0)


    # Fit each of the imputations
    fits <- with(imp, glm(stats::formula(formula_string), family = family))

    # Pool estimates of coefficients
    est <- mice::pool(fits)
    # Grab coefficient summary table
    cf <- summary(est)
    # Extract/name certain columns to make the same output without imputation
    pred_coeff <- cf[2,c(2,3,4,6)]
    colnames(pred_coeff) <- c("Estimate", "Std. Error", "Statistic", "Pr(>|z|)")

  } else {
    # Either imputation is disabled, or imputation is not needed

    # Calculate glm
    mod = stats::glm(data = model_data, stats::formula(formula_string),
                     family = family)

    # Grab predictor from second row of summary output
    cf = summary(mod)$coefficients

    # Get output for predictor
    pred_coeff <- t(as.data.frame(cf[2,]))
    colnames(pred_coeff)[3] <- "Statistic"

  }


  if (verbose) {print(cf)}

  # Produce tidy dataframe of output
  out <- cbind(data.frame(Response = colnames(R)[Rind]),
               Predictor = colnames(P)[Pind],
               pred_coeff,
               data.frame(ADJ_P_VAL = NA, Imputes = n_imputes,
                          Frac_R_Missing = sum(is.na(R[,Rind]))/nrow(R),
                          Frac_P_Missing = sum(is.na(P[,Pind]))/nrow(P),
                          Family = family, Formula = formula_string))
  names(out)[6] <- "P_VAL"
  rownames(out) <- NULL

  return(out)
}
