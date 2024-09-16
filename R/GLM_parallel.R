
#' GLM_parallel
#'
#' @description perform general linear modeling with the equation of two
#' possible forms:
#'
#' Response ~ Parallel_Predictor + Predictors + Confounders
#' Parallel_Response ~ Predictors + Confounders
#'
#' \code{R\[,Rind] ~ P\[,Pind] + Pe\[,1] + Pe\[,...] + C\[, 1] + C\[, ...]}
#'
#' R\[,Rind]: Either a single response variable, or several response variables
#'  where a separate model is fitted for each (parallelized).
#' P\[,Pind]: A series of predictors where a separate model is fitted for each
#'  (parallelized). Note only R or P can be parallelized, but not both.
#' Pe: Extra predictors, to be included in all models.
#' C: confounder variables, to be included in all models.
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
#' @param Pe dataframe of extra predictor variables to be included in all models.
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
#'
#'
#' @importFrom magrittr %>%
#'
GLM_parallel = function(R, Rind = 1, P = NULL, Pind = 1, Pe = NULL, C = NULL,
                        family = "binomial", verbose = FALSE, impute_na = TRUE,
                        db_flag = TRUE) {
  if (db_flag) {save(list = ls(all.names = TRUE), file = "GLM_parallel.RData") }
  # load(file = "GLM_parallel.RData")

  stopifnot("GLM_parallel:error: R must be a dataframe" = is.data.frame(R) ||
              is.matrix(R))
  stopifnot("GLM_parallel:error: P must be a dataframe" = is.data.frame(P) ||
              is.matrix(P) || is.null(P))
  if (length(Rind)!=1 || length(Pind)!=1) {
    stop(sprintf(paste0("GLM_parallel:error: both length() of Rind(==%.0f) and ",
                        "Pind(==%.0f) must be 1."),
                 length(Rind), length(Pind)))
  }

  # If P is null, then set the index to be null also
  if (is.null(P)) Pind = NULL

  # Join all data for model into one dataframe that preserves the variable type
  # for each column, and supports only having some variable types included
  model_data <- R %>% dplyr::select(dplyr::all_of(Rind))
  formula_string <- paste0(colnames(R)[Rind], " ~ ")
  if (!is.null(P)) {
    model_data <- base::cbind(model_data, P %>% dplyr::select(dplyr::all_of(Pind)))
    formula_string <- paste0(formula_string, paste0(colnames(P)[Pind], " + ", collapse = " "))
  }
  if (!is.null(Pe)) {
    model_data <- base::cbind(model_data, Pe)
    formula_string <- paste0(formula_string, paste0(colnames(Pe), c(rep(" + ", max(
      c(ncol(Pe)-1, 0)))), collapse=" "))
  }
  if (!is.null(C)) {
    model_data <- base::cbind(model_data, C)
    formula_string <- paste0(formula_string, paste0(colnames(C), c(rep(" + ", max(
      c(ncol(C)-1, 0)))), collapse=" "))
  }
  formula_string = base::gsub("\\s\\+\\s$", "", formula_string)

  # Check for same value in response
  if ( dim(unique(R[Rind]))[1]==1) {
    save(list = ls(all.names = TRUE), file = "GLM_parallel_error.RData") # load(file = "GLM_parallel.RData")
    stop(sprintf("GLM_parallel: Response variable only has 1 unique value at
                 Rind %i, Pind %i", Rind, Pind))
  }


  # Calculate the number of imputations required
  # Generally, number of imputations is % of missing data
  # Unless no missing data, than imputes is set to 1 to not error package
  # If n_imputes <5, then set to 5, because that is considered the min if there
  # are any missing values
  n_imputes <- ceiling(100 * sum(is.na( model_data[,1:2])) /
                         (2 * nrow(model_data)))
  # If imputation is needed, set to at least 5 (original paper rec)
  if(n_imputes > 0 && n_imputes < 5) {n_imputes = 5}


  # # Start with Response and Predictor
  # formula_string <- paste0(colnames(R)[Rind], " ~ ", colnames(P)[Pind], plus_str(P,is.null)

  # Impute missing data if flag is set and data actually missing
  if (impute_na && n_imputes > 0 ) {

    # Perform multiple imputations of the dataset
    imp <- mice::mice(model_data, print = FALSE, m = n_imputes, maxit = 25, seed = 0)

    # Fit each of the imputations
    fits <- with(imp, glm(stats::formula(formula_string), family = family))

    # Pool estimates of coefficients
    est <- mice::pool(fits)
    # Grab coefficient summary table
    cf <- summary(est)
    rownames(cf) <- cf$term
    cf <- cf[,c(2,3,4,6)]
    colnames(cf) <- c("Estimate", "StdError", "Statistic", "P_VAL")

    aic <- mean(sapply(fits$analyses, function(x) x$aic))

  } else {
    # Either imputation is disabled, or imputation is not needed
    # Calculate glm
    mod = stats::glm(data = model_data, stats::formula(formula_string),
                     family = family)

    # Grab predictor from second row of summary output
    cf = summary(mod)$coefficients

    # If P is not null, remove the extra suffix "TRUE" to its var name
    # Example:  "cg27785526TRUE" ->  "cg27785526"
    rownames(cf)[2] <- stringr::str_replace(rownames(cf)[2], "TRUE$","")

    # cf <- cf[,c(1,2,3, 5)]
    colnames(cf)[2:4] <-c("StdError", "Statistic", "P_VAL")
    aic <- mod$aic

  }

  # Add response column and variable column from rownames
  df_res <-
    cbind(data.frame(Response = colnames(R)[Rind],
                     Variable = rownames(cf)[2:nrow(cf)]),
          cf[2:nrow(cf),,drop = FALSE])
  rownames(df_res) <- NULL

  # Record whether variable is confounder
  # Get number of columns for P[,Pind] and Pe
  n_pred <- max(c(0, length(Pind))) + max(c(0, ncol(Pe)))
  df_res$Confounder <- 0
  if (n_pred > 1) {
      df_res$Confounder[(n_pred +1) : nrow(df_res)] <- rep(1, nrow(df_res) - n_pred)
  }

  # Slot for adjusted p-value, calculated outside of this function
  df_res$ADJ_P_VAL <- NA
  df_res$Family <- family
  df_res$Formula <- formula_string
  df_res$Model_Id <-max(c(Rind, Pind))
  df_res$aic <- aic
  if (verbose) {print(cf)}

  return(df_res)
}
