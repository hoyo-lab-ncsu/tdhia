
#' GLM_parallel
#'
#' @description perform general linear modeling with the equation of the form:
#'
#' Response ~ Predictor + Confounders
#' \code{R[,Rind] ~ P[,Pind] + C[, 1] + C[, 2] + C[, 3] + C[, ...]}
#'
#' One predictor and one response variable is assumed. Either R or P can be
#' parallelized (but only one of them)
#'
#' @param R dataframe where each column is a Response variable.
#' @param Rind column index for R if it has multiple columns (for parallel processing).
#' Default for Rind is 1 for single column dataframe.
#' @param P dataframe where each column is a Predictor variable.
#' @param Pind column index for P if it has multiple columns (for parallel processing).
#' Default for Pind is 1 for single column dataframe.
#' @param C dataframe of confounder variables to be including in model
#' @param family string denoting GLM family
#' @param n_p_adj n adjustment for multiple comparisons
#'
#' @return fit of general linear model, including
#' - Estimate: point estimates of coefficients for each of the predictors
#' - Std. Error: standard error of the estimates
#' - Cumulative two-tailed probability
GLM_parallel = function(R, Rind = 1, P, Pind = 1, C, family = "binomial", n_p_adj = 1) {
  save(list = ls(all.names = TRUE), file = "GLM_parallel.RData")
  # load(file = "GLM_parallel.RData")

  stopifnot("GLM_parallel:error: R must be a dataframe" = is.data.frame(R) || is.matrix(R))
  stopifnot("GLM_parallel:error: P must be a dataframe" = is.data.frame(P) || is.matrix(P))
  if (length(Rind)!=1 || length(Pind)!=1) {
    stop(sprintf("GLM_parallel:error: both length() of Rind(==%.0f) and Pind(==%.0f) must be 1.",
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
  mod = stats::glm(stats::formula(formula_string), family = family)

  # Grab predictor from second row of summary output
  cf = summary(mod)$coefficients

  # Produce tidy dataframe of output
  out <- cbind(data.frame(Response = colnames(R)[Rind]),
               Predictor = colnames(P)[Pind],
               t(as.data.frame(cf[2,])),
               data.frame(ADJ_P_VAL = stats::p.adjust(cf[2,4], method = "fdr", n = n_p_adj),
                 Formula = formula_string, Family = family))
  names(out)[6] <- "P_VAL"

  rownames(out) <- NULL

  return(out)
}
