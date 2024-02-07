

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
#' values. The code handles either case. Note that R and P cannot both be
#' multicolumn.
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
#' @param impute_na boolean flag, when TRUE missing values in R or P are
#' recovered with multiple imputations with the mice package.
#' @param n.cores number of cores used for parallel computation. Default is NULL
#' and set to max cores detected minus 1.
#' @param db_flag boolean when true export workspace to disk for debugging.
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %dopar%
#' @importFrom rlang .data
#'
#' @return results with fitted coefficients from the glm, sorted by p-value
study_imprint <- function (R, P, C, family, n_p_adj = max(c(ncol(R), ncol(P))),
                           max_p_val = 0.05, impute_na = TRUE, n.cores = NULL,
                           db_flag = FALSE) {

  if (db_flag) {save(list = ls(all.names = TRUE), file = "study_imprint.RData")}
  # load(file = "study_imprint.RData")

  cat("Analyzing imprintome study...\n")

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

  # Calculate missing Values
  fract_r_na <- sum(is.na(R))/(nrow(R)*ncol(R))
  fract_p_na <- sum(is.na(P))/(nrow(P)*ncol(P))
  cat(sprintf("  R: %.0f%% of values are missing (NA values).\n", fract_r_na*100 ))
  cat(sprintf("  P: %.0f%% of values are missing (NA values).\n", fract_p_na*100 ))
  n_imput <- mean(c(fract_r_na, fract_p_na))

  #  Get indexes for response and predictor variables
  Rind <- 1:ncol(R)
  Pind <- 1:ncol(P)

  # Print short label for model (excluding confounders)
  cat(sprintf("First model label: %s ~ %s + confounders\n", colnames(R)[1],
              colnames(P)[1]))
  cat("Example output of first model:\n")

  # Define parallel processing function for GLM
  #-----------------------------------------------------------------------------
  # Run function on first input in verbose to show model for debugging
  if (length(Rind) > 1) {
    # For each response variable
    foreach_fun <- function (x) GLM_parallel(R = R, Rind = x, P = P, Pind = Pind,
                                             C = C, family = family,
                                             impute_na = impute_na,
                                             db_flag = db_flag)
    test_run <- GLM_parallel(R = R, Rind = 1, P = P, Pind = Pind,
                             C = C, family = family, verbose = TRUE)
  } else if (length(Pind) > 1) {
    # For each predictor variable
    foreach_fun <-  function (x) GLM_parallel(R = R, Rind = Rind, P = P, Pind = x,
                                              C = C, family = family,
                                              impute_na = impute_na,
                                              db_flag = db_flag)
    test_run <- GLM_parallel(R = R, Rind = Rind, P = P, Pind = 1,
                             C = C, family = family, verbose = TRUE)
  } else {
    foreach_fun <- function (x) GLM_parallel(R = R, Rind = Rind, P = P, Pind = Pind,
                                             C = C, family = family,
                                             impute_na = impute_na,
                                             db_flag = db_flag)
    test_run <- GLM_parallel(R = R, Rind = Rind, P = P, Pind = Pind,
                             C = C, family = family, verbose = TRUE)
  }
  cat(sprintf("Formula: %s\n", test_run$Formula))


  # Identify parallel index processing
  if (length(Rind)>1) {par_ind <- Rind
  }  else if (length(Pind)>1) {par_ind <- Pind
  }  else {par_ind <- 1}

  cat(sprintf(">>  Fitting all %.0f models...\n", length(par_ind)))
  # Fit models to data in parallel
  start <- Sys.time()

  # Initialize parallel computing
  if (is.null(n.cores)) { n.cores = parallel::detectCores() - 1}
  cl <- snow::makeCluster(n.cores)
  doSNOW::registerDoSNOW(cl)

  # Add progress bar
  pb <- utils::txtProgressBar(max = max(par_ind), style = 3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  # Execute parallel processing
  df_fits <- foreach::foreach(x = par_ind, .combine = rbind,
                              .export = "GLM_parallel",.options.snow = opts
                              ) %dopar%
    {
      foreach_fun(x)
    }

  # End processing time
  finish <- Sys.time()
  cat(" Finished.\n")
  # cat(sprintf(">Processing time: %f %s\n", finish - start, units(finish - start)))

  # Sort by p value
  df_fits_sorted <- df_fits %>% dplyr::arrange(.data$P_VAL)
  df_fits_sorted$P_VAL[is.na(df_fits_sorted$P_VAL)] <- 1
  # Recalculate P adjust
  df_fits_sorted$ADJ_P_VAL <- custom_p.adjust(df_fits_sorted$P_VAL, method = "fdr",
                                              n = n_p_adj)


  # Report results from GLM
  cat(sprintf(">>  %.0f cpg sites have p_val < %.2f\n",
              sum(df_fits_sorted$P_VAL < max_p_val), max_p_val))
  cat(sprintf(">>  %.0f cpg sites have adj_p_val < %.2f\n",
              sum(df_fits_sorted$ADJ_P_VAL < max_p_val), max_p_val))
  print(df_fits_sorted[df_fits_sorted$ADJ_P_VAL < max_p_val, c(1:7,11)])
  cat("\n\n")


  # TODO: fill in what this does
  lambda <-
    stats::median(stats::qchisq( as.numeric( as.character(df_fits_sorted$P_VAL)),
                                 df = 1, lower.tail = F),
                  na.rm = T) / stats::qchisq( 0.5, 1)
  # cat(sprintf("> Lambda: %f\n", lambda))


  return(df_fits_sorted)
}







