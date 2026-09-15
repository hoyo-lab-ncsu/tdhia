
#' Fit One Model in a Series of GLM Association Tests
#'
#' Selects a response and primary predictor, adds common covariates, and
#' fits one generalized linear model. This worker can be called repeatedly
#' or in parallel by an outer analysis function.
#'
#' @param R Data frame or matrix of response variables with samples in rows.
#' @param Rind Single column index selecting the response in R.
#' @param P Optional data frame or matrix of primary predictors. Rows must
#'   already correspond to the same samples and order as R.
#' @param Pind Single column index selecting the primary predictor in P.
#'   Ignored when P is NULL.
#' @param Pe Optional data frame of additional predictors included together.
#'   Rows must already be aligned with R.
#' @param family GLM family specification passed to .fit_model.
#' @param verbose Logical; print the fitted coefficient table.
#' @param impute_na Logical; enable multiple imputation through .fit_model.
#' @param db_flag Logical; save the initial environment to GLM_parallel.RData
#'   in the working directory. Failure snapshots can be written even if FALSE.
#' @param .fit_model Fitting function with arguments model_data,
#'   formula_string, family, impute_na, and n_imputes, in that order.
#'   Must return a list containing cf (a coefficient table including an
#'   intercept row) and aic. Defaults to fit_model().
#'
#' @details
#' The model uses the selected R column as response, the selected P column
#' when supplied, and all Pe columns as predictors. Sample alignment is the
#' caller's responsibility. The number of imputations is the ceiling of the
#' percentage of missing model-data cells, with a minimum of five when any
#' values are missing.
#'
#' Errors raised by .fit_model are caught. A failed fit yields placeholder
#' coefficients with NA estimates and p-values of 1, and writes a
#' GLM_parallel_error_pid-<process ID>.RData snapshot. A constant response
#' detected before fitting writes GLM_parallel_error.RData and stops.
#'
#' @return A data frame with one row per non-intercept coefficient:
#'   - Response and Variable: response and coefficient names.
#'   - Estimate, StdError, Statistic, P_VAL: coefficient summary values.
#'   - Confounder: legacy indicator initialized to 0; the current single-index
#'     implementation does not mark additional covariates as confounders.
#'   - ADJ_P_VAL: NA placeholder for adjustment by the caller.
#'   - Family and Formula: family specification and constructed formula.
#'   - Model_Id: maximum of the selected response and predictor indices.
#'   - aic: model AIC, or mean AIC for imputed fits.
#' @seealso [fit_model()], [imprintome_glm()]
#' @importFrom magrittr %>%
#' @export
GLM_parallel = function(R, Rind = 1, P = NULL, Pind = 1, Pe = NULL,
                        family = "binomial", verbose = FALSE, impute_na = TRUE,
                        db_flag = TRUE, .fit_model = fit_model) {
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
  # n_imputes <- ceiling(100 * sum(is.na( model_data[,1:2])) /
  #                        (2 * nrow(model_data)))
  n_imputes <- ceiling(sum(is.na(model_data)) / (nrow(model_data)*ncol(model_data)) *100)
  # If imputation is needed, set to at least 5 (original paper rec)
  if(n_imputes > 0 && n_imputes < 5) {n_imputes = 5}

  
  # Fit model and return fits and AIC, if fail save output and empty results
  coefs <- stringr::str_extract_all(formula_string, stringr::regex("[^\\s-\\+~]+"))[[1]][-1]
  cf = data.frame(Estimate = rep(NA, length(coefs)+1), StdError = NA, 
                  Statistic = NA, P_VAL = 1, row.names = c("(Intercept)",coefs))
  aic <- NA; model_out <-  NULL
  tryCatch({
    # fdsdf(fdsfsd)
    model_out <- .fit_model(model_data, formula_string, family, impute_na, n_imputes)
    cf = model_out$cf
    aic = model_out$aic
 },  error = function(e) {})
 if (is.null(model_out)) {
   save(list = ls(all.names = TRUE), file = paste0("GLM_parallel_error_pid-", Sys.getpid(),".RData"))
 }
  
  # Add response column and variable column from rownames
  df_res <-
    cbind(data.frame(Response = colnames(R)[Rind],
                     Variable = rownames(cf)[2:nrow(cf)]),
          cf[2:nrow(cf),,drop = FALSE])
  rownames(df_res) <- NULL

  # Record whether variable is confounder
  # Get number of columns for P[,Pind] and Pe
  n_pred <- max(c(0, length(Pind)))
  df_res$Confounder <- 0
  if (n_pred > 1) {
      df_res$Confounder[(n_pred) : nrow(df_res)] <- rep(1, nrow(df_res) - n_pred)
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


#' Fit a GLM with Optional Multiple Imputation
#'
#' Fits a generalized linear model directly or fits and pools models across
#' MICE imputations, returning a common coefficient-table format.
#'
#' @param model_data Data frame containing all variables in formula_string.
#' @param formula_string Character string specifying a model formula.
#' @param family GLM family specification accepted by glm(), usually
#'   "gaussian" or "binomial".
#' @param impute_na Logical; use multiple imputation when TRUE and
#'   n_imputes is greater than zero.
#' @param n_imputes Number of imputations passed as m to mice::mice().
#'   Zero bypasses imputation even when impute_na is TRUE.
#'
#' @details
#' The imputation branch uses maxit = 20 and seed = 0, then pools coefficient
#' estimates with mice::pool(). The direct branch uses the usual glm()
#' missing-data handling. Both branches remove a trailing TRUE suffix from
#' the first non-intercept coefficient's name.
#'
#' @return A named list with cf and aic. The cf table includes the intercept
#'   and has coefficient names as row names, with columns Estimate, StdError,
#'   Statistic, and P_VAL. The aic value is the fitted model's AIC or the mean
#'   AIC across imputed-data models. Statistic and P_VAL come from the direct
#'   GLM summary or pooled coefficient summary, respectively.
#' @seealso [GLM_parallel()]
#' @export
fit_model <- function(model_data, formula_string, family, impute_na, n_imputes) {
  
  
  # Impute missing data if flag is set and data actually missing
  if (impute_na && n_imputes > 0 ) {
    
    # Perform multiple imputations of the dataset
    imp <- mice::mice(model_data, print = FALSE, m = n_imputes, maxit = 20, seed = 0)
    
    # Fit each of the imputations
    fits <- with(imp, glm(stats::formula(formula_string), family = family))
    # Pool estimates of coefficients
    est <- mice::pool(fits)
    # Grab coefficient summary table
    cf <- summary(est)
    rownames(cf) <- cf$term
    cf <- cf[,c(2,3,4,6)]
    rownames(cf)[2] <- stringr::str_replace(rownames(cf)[2], "TRUE$","")
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
    colnames(cf)[2:4] <-c("StdError", "Statistic", "P_VAL")
    aic <- mod$aic
    
  }
  return(list(cf = cf, aic = aic))
}