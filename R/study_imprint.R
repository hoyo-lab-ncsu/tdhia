

#' study_imprint
#'
#' @description fits a series of general linear models with parallel processing
#' using response variable(s) R, predictor variable(s) R and confouding
#' variable(s) C. The model formula takes on the form:
#'
#' If R has many variables/ columns and P has one, models are parallelized over R:
#'
#' R\[,i] ~ P\[,1] + C\[,1] + C\[,2] + C\[,3] ...
#'
#' If P has many variables/ columns and R has one, models are parallelized over P:
#'
#' R\[,1] ~ P\[,i] + C\[,1] + C\[,2] + C\[,3] ...
#'
#' This designed that either R or P are the CpG beta values or the ICR beta
#' values. The code handles other cases.
#'
#' @param R response variable(s): data frame that is patients (rows) x variables
#'  (columns) could either be beta values for cpg or ICR sites, or
#'  study metadata.
#' @param P predictor variable(s) data frame that is patients (rows) x variables
#'  (columns), could either be beta values for cpg or ICR sites, or
#'  study metadata.
#' @param C predictor variable(s) data frame that is patients (rows) x variables
#'  (columns).
#' @param family string for family model argument for glm(). (ex. 'binomial', 'gaussian')
#' @param n_p_adj number of comparisons to use for p-value adjustment. Sometimes
#'  is this is different than the number of actual comparisons (could be more or
#'  less). Defaults to max columns to response and predictor variables
#' @param max_p_val max p-value threshold, prints number of entries below this
#'  threshold, but does not alter computation.
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %dopar%
#' @importFrom rlang .data
#'
#' @return results with fitted coefficients from the glm, sorted by p-value
study_imprint <- function (R, P, C, family, n_p_adj = max(c(ncol(R), ncol(P))),
                                      max_p_val = 0.05) {

  save(list = ls(all.names = TRUE), file = "analyze_imprintome_study.RData")
  # load(file = "analyze_imprintome_study.RData")

  # Reorder rows in P and C to match R, keep both as dataframes even if 1 column
  P <- P %>%
    tibble::rownames_to_column(var = "row_names") %>%
    dplyr::arrange(factor(.data$row_names, levels = rownames(R))) %>%
    tibble::column_to_rownames(var = "row_names")
  if( !all(rownames(R)==rownames(P))) {
    stop("glm: rownames of R and P dataframes do not match.")}

  C <- C %>%
    tibble::rownames_to_column(var = "row_names") %>%
    dplyr::arrange(factor(.data$row_names, levels = rownames(R))) %>%
    tibble::column_to_rownames(var = "row_names")
  if( !all(rownames(R)==rownames(C))) {
    stop("glm: rownames of R and C dataframes do not match.")}


  #  Get indexes for response and predictor variables
  Rind <- 1:ncol(R)
  Pind <- 1:ncol(P)

  # Print short label for model (excluding confounders)
  cat(sprintf("Model Label: %s ~ %s\n", colnames(R)[1], colnames(P)[1]))
  cat("Output of first model:\n")

  # Define parallel processing function for GLM
  # Run function on first input in verbose to show model for debugging
  if (length(Rind) > 1) {
    # For each response variable
    foreach_fun <- function (x) GLM_parallel(R = R, Rind = x, P = P, Pind = Pind,
                                             C = C, family = family)
    test_run <- GLM_parallel(R = R, Rind = 1, P = P, Pind = Pind,
                 C = C, family = family, verbose = TRUE)
  } else if (length(Pind) > 1) {
    # For each predictor variable
    foreach_fun <-  function (x) GLM_parallel(R = R, Rind = Rind, P = P, Pind = x,
                                              C = C, family = family)
    test_run <- GLM_parallel(R = R, Rind = Rind, P = P, Pind = 1,
                 C = C, family = family, verbose = TRUE)
  } else {
    foreach_fun <- function (x) GLM_parallel(R = R, Rind = Rind, P = P, Pind = Pind,
                                             C = C, family = family, )
    test_run <- GLM_parallel(R = R, Rind = Rind, P = P, Pind = Pind,
                 C = C, family = family, verbose = TRUE)
  }
  cat(sprintf("Formula: %s\n", test_run$Formula))


  # Initialize parallel computing
  cl <- parallel::makePSOCKcluster(parallel::detectCores() - 1)
  doParallel::registerDoParallel(cl)

  # Identify parallel index processing
  if (length(Rind)>1) {par_ind <- Rind
  }  else if (length(Pind)>1) {par_ind <- Pind
  }  else {par_ind <- 1}

  # Fit mdoels to data in parallel
  start <- Sys.time()
  cat(sprintf(">>  Fitting all %.0f models...", length(par_ind)))
  df_fits <- foreach::foreach(x = par_ind, .combine = rbind,
                              .export = "GLM_parallel") %dopar% {
                                foreach_fun(x)
                              }
  # End processing time
  finish <- Sys.time()
  cat(" Finished.\n")
  # cat(sprintf(">Processing time: %f %s\n", finish - start, units(finish - start)))

  # Sort by p value
  df_fits_sorted <- df_fits %>% dplyr::arrange(.data$P_VAL)
  # Recalculate P adjust
  df_fits_sorted$ADJ_P_VAL <- custom_p.adjust(df_fits_sorted$P_VAL, method = "fdr",
                                              n = n_p_adj)

  # Report results from GLM
  cat(sprintf(">>  %.0f cpg sites have p_val < %.2f\n",
              sum(df_fits_sorted$P_VAL < max_p_val), max_p_val))
  cat(sprintf(">>  %.0f cpg sites have adj_p_val < %.2f\n",
              sum(df_fits_sorted$ADJ_P_VAL < max_p_val), max_p_val))
  print(df_fits_sorted[df_fits_sorted$ADJ_P_VAL < max_p_val, c(1:7,9)])
  cat("\n\n")


  # TODO: fill in what this does
  lambda <-
    stats::median(stats::qchisq( as.numeric( as.character(df_fits_sorted$P_VAL)),
                                 df = 1, lower.tail = F),
                  na.rm = T) / stats::qchisq( 0.5, 1)
  # cat(sprintf("> Lambda: %f\n", lambda))


  return(df_fits_sorted)
}







