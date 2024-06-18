

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
#' @param P parallelized predictor variable(s) data frame that is patients
#'  (rows) x variables (columns), could either be beta values for cpg or ICR
#'  sites. For study metadata, use Pe input argument.
#' @param Pe extra predictor variables that are not parallelized (incuded in all
#'  fitted models).
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
#' @param rm.na.R flag to remove NA values in R prior to imputation.
#' @param rm.na.P flag to remove NA values in P prior to imputation.
#' @param rm.na.Pe flag to remove NA values in Pe prior to imputation.
#' @param rm.na.C flag to remove NA values in C prior to imputation.
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %dopar%
#' @importFrom rlang .data
#'
#' @return results with fitted coefficients from the glm, sorted by p-value
study_imprint <- function (R, P, Pe, C, family, n_p_adj = max(c(ncol(R), ncol(P))),
                           max_p_val = 0.05, impute_na = TRUE, n.cores = NULL,
                           db_flag = FALSE, rm.na.R = FALSE, rm.na.P = FALSE,
                           rm.na.Pe = FALSE, rm.na.C = FALSE, rm.na.all = FALSE,
                           print_confounders = FALSE, verbose = TRUE,
                           icr_mapping = NULL) {
  if (db_flag) {save(list = ls(all.names = TRUE), file = "study_imprint.RData")}
  # load(file = "study_imprint.RData")

  if (rm.na.all) rm.na.R <- rm.na.P <- rm.na.Pe <- rm.na.Pe <- rm.na.C <- TRUE
  # This expression prevents devtools from issuing a NOTE warning
  # x is defined within some local functions below
  x <- NULL
  verbosecat = function (x) if (verbose) cat(x)

  verbosecat("Analyzing imprintome study...\n")

  # Reorder rows between P,Pe,C to match R
  # ----------------------------------------------------------------------------
  if (!is.null(P)) {
  P <- P %>%
    tibble::rownames_to_column(var = "row_names") %>%
    dplyr::arrange(factor(.data$row_names, levels = rownames(R))) %>%
    tibble::column_to_rownames(var = "row_names")
  if( !all(rownames(R)==rownames(P))) {
    stop("glm: rownames of R and P dataframes do not match.")}
  }

  if (!is.null(Pe)) {
    Pe <- Pe %>%
      tibble::rownames_to_column(var = "row_names") %>%
      dplyr::arrange(factor(.data$row_names, levels = rownames(R))) %>%
      tibble::column_to_rownames(var = "row_names")
    if( !all(rownames(R)==rownames(Pe))) {
      stop("glm: rownames of R and Pe dataframes do not match.")}
  }

  if (!is.null(C)) {
    C <- C %>%
      tibble::rownames_to_column(var = "row_names") %>%
      dplyr::arrange(factor(.data$row_names, levels = rownames(R))) %>%
      tibble::column_to_rownames(var = "row_names")
    if( !all(rownames(R)==rownames(C))) {
      stop("glm: rownames of R and C dataframes do not match.")}
  }

  # Calculate missing Values
  fract_r_na <- sum(is.na(R))/(nrow(R)*ncol(R))
  fract_p_na <- sum(is.na(P))/(nrow(P)*ncol(P))
  verbosecat(sprintf("  R: %.0f%% of values are missing (NA values).\n", fract_r_na*100 ))
  verbosecat(sprintf("  P: %.0f%% of values are missing (NA values).\n", fract_p_na*100 ))
  n_imput <- mean(c(fract_r_na, fract_p_na))

  # Remove rows in dataset with NA values (if requested)
  # NAs in parallel variable (for cpgs or icrs) can be imputed in model fitting
  #-----------------------------------------------------------------------------
  verbosecat("Checking for NAs and Removing if specified.\n")

  is.R.na = rep(FALSE, nrow(R))
  if (!is.null(R)) {
    temp <- rowSums(is.na(as.matrix(R))) > 0
    verbosecat(sprintf("   rm.na, R: %.0f/ %.0f rows have 1+ NAs...", sum(temp),
                length(is.R.na)))
    if (rm.na.R) { is.R.na = temp
      verbosecat(" Marked for removal.\n")
    } else {verbosecat("Keeping them.\n")}
  }

  is.P.na = rep(FALSE, nrow(R))
  if (!is.null(P)) {
    temp =  rowSums(is.na(as.matrix(P))) > 0
    verbosecat(sprintf("   rm.na, P: %.0f/ %.0f rows have 1+ NAs...", sum(temp),
                length(is.P.na)))
    if (rm.na.P) { is.P.na = temp
    verbosecat(" Marked for removal.\n")
    } else {verbosecat("Keeping them.\n")}
  }

  is.Pe.na = rep(FALSE, nrow(R))
  if (!is.null(Pe)) {
    temp = rowSums(is.na(as.matrix(Pe))) > 0
    verbosecat(sprintf("   rm.na, Pe: %.0f/ %.0f rows have 1+ NAs...", sum(temp),
                length(is.Pe.na)))
    if (rm.na.Pe) { is.Pe.na = temp
    verbosecat(" Marked for removal.\n")
    } else {verbosecat("Keeping them.\n")}
  }

  is.C.na = rep(FALSE, nrow(R))
  if (!is.null(C)) {
    temp = rowSums(is.na(as.matrix(C))) > 0
    verbosecat(sprintf("   rm.na, C: %.0f/ %.0f rows have 1+ NAs...", sum(temp),
                length(is.C.na)))
    if (rm.na.C) { is.C.na = temp
    verbosecat(" Marked for removal.\n")
    } else {verbosecat("Keeping them.\n")}
  }


  # Remove all rows marked for removal
  rm.na.flags <- is.R.na | is.P.na | is.Pe.na | is.C.na
  verbosecat(sprintf(">> rm.na: Removing %.0f rows total (before imputation)...\n",
              sum(rm.na.flags)))
  verbosecat(sprintf(">> rm.na:%.0f rows now remain.\n",
              sum(!rm.na.flags)))

  R  <- R[!rm.na.flags, , drop = FALSE]
  P  <- P[!rm.na.flags, , drop = FALSE]
  Pe <- Pe[!rm.na.flags, , drop = FALSE]
  C  <- C[!rm.na.flags, , drop = FALSE]


  #  Get indexes for response and predictor variables
  #-----------------------------------------------------------------------------
  Rind <- 1: max(is.null(ncol(R)),ncol(R))
  Pind <- 1: max(is.null(ncol(P)),ncol(P))

  # Print short label for model (excluding confounders)
  verbosecat(sprintf("First model label: %s ~ %s ...\n", colnames(R)[1],
              colnames(P)[1]))
  verbosecat("Example output of first model:\n")

  # Define parallel processing function for GLM
  #-----------------------------------------------------------------------------
  # Run function on first input in verbose to show model for debugging
  if (length(Rind) > 1) {
    # For each response variable
    foreach_fun <- function (x) GLM_parallel(R = R, Rind = x, P = P, Pind = Pind,
                                             Pe = Pe, C = C, family = family,
                                             impute_na = impute_na,
                                             db_flag = db_flag, verbose = verbose)
    test_run <- GLM_parallel(R = R, Rind = 1, P = P, Pind = Pind,
                             Pe = Pe, C = C, family = family,
                             impute_na = impute_na,
                             db_flag = db_flag, verbose = verbose)
  } else if (length(Pind) > 1) {
    # For each predictor variable
    foreach_fun <-  function (x) GLM_parallel(R = R, Rind = Rind, P = P, Pind = x,
                                              Pe = Pe, C = C, family = family,
                                              impute_na = impute_na,
                                              db_flag = db_flag, verbose = verbose)
    test_run <- GLM_parallel(R = R, Rind = Rind, P = P, Pind = 1,
                             Pe = Pe, C = C, family = family,
                             impute_na = impute_na,
                             db_flag = db_flag, verbose = verbose)
  } else {
    foreach_fun <- function (x) GLM_parallel(R = R, Rind = Rind, P = P, Pind = Pind,
                                             Pe = Pe, C = C, family = family,
                                             impute_na = impute_na,
                                             db_flag = db_flag, verbose = verbose)
    test_run <- GLM_parallel(R = R, Rind = 1, P = P, Pind = 1,
                             Pe = Pe, C = C, family = family,
                             impute_na = impute_na,
                             db_flag = db_flag, verbose = verbose)
  }
  verbosecat(sprintf("Formula: %s\n", test_run$Formula))


  # Identify parallel index processing
  if (length(Rind)>1) {par_ind <- Rind
  }  else if (length(Pind)>1) {par_ind <- Pind
  }  else {par_ind <- 1}

  verbosecat(sprintf(">>  Fitting all %.0f models...\n", length(par_ind)))
  # Fit models to data in parallel
  start <- Sys.time()

  # Initialize parallel computing
  if (is.null(n.cores)) { n.cores = parallel::detectCores() - 1}
  cl <- snow::makeCluster(n.cores)
  doSNOW::registerDoSNOW(cl)

  # Add progress bar
  # if (verbose) {
    pb <- utils::txtProgressBar(max = max(par_ind), style = 3)
    opts <- list(progress =  function(n) utils::setTxtProgressBar(pb, n))
  # } else { opts = NULL }
  # Execute parallel processing
  df_fits <- foreach::foreach(x = par_ind, .combine = rbind,
                              .export = "GLM_parallel",
                              .packages = c("magrittr"),
                              .options.snow = opts
  ) %dopar% {
    foreach_fun(x)
  }
  # End processing time
  cat("\n")
  finish <- Sys.time()
  verbosecat(" Finished.\n")


  # Separate results for each variable of model
  verbosecat("Separating results for each variable used in model...\n")
  model_vars <-   df_fits$Variable[df_fits$Model_Id==1]
  if (!is.null(P) ) {model_vars <- model_vars[2:length(model_vars)]}
  model_vars <- model_vars[!is.na(model_vars)]

  # Go through each model variable and extract results
  dfs_sep <- list()
  for (n in seq_along(model_vars)) {
    dfs_sep[[n]] <- df_fits[df_fits$Variable == model_vars[n],]
  }
  names(dfs_sep) <- model_vars
  # If cpg or ICR sites are a predictor, name the dataframe "site_id
  if (!is.null(P)) {
    df_temp <- list(df_fits[df_fits$Variable %in% colnames(P),])
    names(df_temp) <- "imp_site"
    dfs_sep <- c(df_temp, dfs_sep)
  }

  # Set NA p-values to a max value of 1 (for sorting)
  na_fun <- function(df) { df$P_VAL[is.na(df$P_VAL)] <- 1; return(df)}
  # Sort each variable by p-value
  sort_fun <- function(df) {df=df[order(df$P_VAL, decreasing = FALSE),]; return(df)}


  # Export for debugging
  if (db_flag) {save(list = ls(all.names = TRUE), file = "study_imprint2.RData")}
  # load(file = "study_imprint2.RData")


  # Calculate adjusted p_value
  verbosecat("Adjusting p-values...\n")
  adj_p_val <- function(df) {
    df$ADJ_P_VAL <- custom_p.adjust(df$P_VAL, method = "fdr", n = n_p_adj);
    return(df)
  }
  dfs_sorted <- lapply(dfs_sep, na_fun)
  dfs_sorted <- lapply(dfs_sorted, sort_fun)
  dfs_corr   <- lapply(dfs_sorted, adj_p_val)

  dfs_corr$example_formula <- test_run$Formula[1]

  # Print out results of analysis
  summarize_study(dfs_corr, varname = "imp_site",max_p_val, print_sites = FALSE,
                  print_confounders = FALSE)


  return(dfs_corr)
}



# TODO: fill in what this does
# lambda <-
#   stats::median(stats::qchisq( as.numeric( as.character(df_fits_sorted$P_VAL)),
#                                df = 1, lower.tail = F),
#                 na.rm = T) / stats::qchisq( 0.5, 1)
# verbosecat(sprintf("> Lambda: %f\n", lambda))



