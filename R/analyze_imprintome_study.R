

#' analyze_imprintome_study
#'
#' @description performs a statistical analysis with CpG sites whose methylation
#' state are correlated with a response variable while also considering several
#' covariates. CpG methylation status/ beta values are used as predictors, from
#' the study metadata a response variable is chosen along with any co-factors. The
#' presence of an association is tested for each of the predictor variables to
#' the response variable.
#'
#' @param model_params a named list that defines parameter for the model to be
#' fitted to the data, with the following fields:
#'  - R: response variable(s), data frame that is patients (rows) x variables
#'  (columns) could either be beta values for cpg or ICR sites, or
#'  study metadata.
#'  - Rind: single value or array that represents column index for which variables
#'  in R to use in the model, for the case of multiple columns, will be parallelized
#'  in separate models. Only Rind or Pind can have more than one value, but not
#'  both (defaults to 1).
#'  - P: predictor variable(s) data frame that is patients (rows) x variables
#'  (columns), could either be beta values for cpg or ICR sites, or
#'  study metadata.
#'  - Pind: single value or array that represents column index for which variables
#'  in P to use in the model, for the case of multiple columns, will be parallelized
#'  in separate models. Only Rind or Pind can have more than one value, but not
#'  both (defaults to 1).
#'  - C: confounder variables data frame that is patients (rows) x variables
#'  (columns).
#'  - family: string for family model argument for glm(). (ex. 'binomial', 'gaussian')
#'  @param quantile_norm boolean flag, when true the beta matrix is quantile normalized
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %dopar%
#' @importFrom rlang .data
#'
#' @return sorted results from statistical analysis
#'
analyze_imprintome_study <- function (model_params) {

  save(list = ls(all.names = TRUE), file = "analyze_imprintome_study.RData")
  # load(file = "analyze_imprintome_study.RData")

  # Export model parameters while keeping both as dataframes even if 1 column
  # Response variable(s)
  R <- model_params[["R"]]
  # Predictor variable(s)
  P <- model_params[["P"]]
  # Confounder Variable(s)
  C <- model_params[["C"]]
  family <- model_params$family

  # Reorder rows in P and C to match R, keep both as dataframes even if 1 column
  P <- P %>%
    tibble::rownames_to_column(var = "row_names") %>%
    dplyr::arrange(factor(.data$row_names, levels = rownames(R))) %>%
    tibble::column_to_rownames(var = "row_names")
  C <- C %>%
    tibble::rownames_to_column(var = "row_names") %>%
    dplyr::arrange(factor(.data$row_names, levels = rownames(R))) %>%
    tibble::column_to_rownames(var = "row_names")

  if( !all(rownames(R)==rownames(P)))
  {stop("glm: rownames of R and P dataframes do not match.")}
  if( !all(rownames(R)==rownames(C)))
  {stop("glm: rownames of R and C dataframes do not match.")}


  #  Get indexes for response and predictor variables
  Rind <- 1:ncol(R)
  Pind <- 1:ncol(P)

  # Helper function for parallel processing of fitting models, either parallelize
  # over columns in R or P
  if (length(Rind) > 1) {
    # For each response variable
    foreach_fun <- function (x) GLM_parallel(R = R, Rind = x,
                                             P = P, Pind = Pind,
                                             C = C, family = family)
  } else if (length(Pind) > 1) {
    # For each predictor variable
    foreach_fun <-  function (x) GLM_parallel(R = R, Rind = Rind,
                                              P = P, Pind = x,
                                              C = C, family = family)
  } else {
    foreach_fun <- function (x) GLM_parallel(R = R, Rind = Rind,
                                             P = P, Pind = Pind,
                                             C = C, family = family)
  }
  cat("Attempting single run of model for debugging...\n")
  foreach_fun(1)

  #Initialize parallel computing
  cl <- parallel::makePSOCKcluster(parallel::detectCores() - 1)
  doParallel::registerDoParallel(cl)

  # Identify parallel index processing
  if (length(Rind)>1) {par_ind <- Rind
  }  else if (length(Pind)>1) {par_ind <- Pind
  }  else {par_ind <- 1}

  # Fit mdoels to data in parallel
  start <- Sys.time()
  cat("Fitting model...\n")
  df_fits <- foreach::foreach(x = par_ind, .combine = rbind,
                              .export = "GLM_parallel") %dopar% {
                                foreach_fun(x)
                              }
  # Rename Columns
  colnames(df_fits) <-  c("Response", "Predictor", "EST", "SE", "Z", "P_VAL", "Formula", "Family")
  # Adjust p-values for FDR correction
  df_fits$ADJ_P_VAL <- stats::p.adjust(df_fits$P_VAL, method = "fdr")
  # Reorder columns for readability
  df_fits <- df_fits[,c(1:6,9, 7, 8)]

  # End processing time
  finish <- Sys.time()
  cat(sprintf("Processing time: %f %s\n", finish - start, units(finish - start)))

  # Sort by p value
  df_fits_sorted <- df_fits %>% dplyr::arrange(.data$P_VAL)

  # Report results from GLM
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







