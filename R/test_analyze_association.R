

parse_lm_formula <- function(formula) {
  # Accept either a character string or a formula
  if (is.character(formula)) {
    formula <- stats::as.formula(formula)
  }
  
  response <- all.vars(formula[[2]])
  
  terms_obj <- stats::terms(formula)
  
  # Labels of all RHS terms
  rhs_terms <- attr(terms_obj, "term.labels")
  
  primary_predictor <- if (length(rhs_terms) >= 1) rhs_terms[1] else NULL
  
  covariates <- if (length(rhs_terms) > 1) rhs_terms[-1] else character(0)
  
  list(
    response = response,
    primary_predictor = primary_predictor,
    covariates = covariates
  )
}



#' imprintome_glm
#'
#' @description fits a series of general linear models with parallel processing
#'  across each CpG or ICR site of the form
#'  
#'  R ~ P + Pe
#'  Where R is the response variable, P is the primary predictor, and Pe are the
#'  covariate variables. The beta matrix (either CpG or ICR level) must be 
#'  specified with the term "beta" in the model string.
#'
#' @param model_str model string for glm. Should use column names in study_data, 
#' and the term "beta" to reference the beta methylation data. Beta must either
#'  be a response or primary predictor.
#' @param study_data dataframe of study variables (each row is a sample). 
#' Rownames should be patient labels that match columns in betas.
#' @param betas beta matrix, each row is imprintome site (cpg or ICR), each 
#' column is a sample (patient). Column names should be patient labels that match 
#' rownames in study_data.
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
#' @param rm.na.all flag to remove all NA values across R, P, and Pe
#' @param verbose prints additional info to console
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %dopar%
#' @importFrom rlang .data
#'
#' @return results with fitted coefficients from the glm, sorted by p-value
#' @export
imprintome_glm <- function (
    model_str, study_data, betas, family, n_p_adj = NULL, max_p_val = 0.05, impute_na = TRUE, n.cores = NULL,
    db_flag = FALSE, rm.na.R = FALSE, rm.na.P = FALSE, m_value_transform = T,
    rm.na.Pe = FALSE, rm.na.all = FALSE, verbose = TRUE) {
  
  
  if (db_flag) {save(list = ls(all.names = TRUE), file = "imprintome_glm.RData")}
  # load(file = "imprintome_glm.RData")
  
  # This expression prevents devtools from issuing a NOTE warning
  # x is defined within some local functions below
  x <- NULL
  verbosecat = function (x) if (verbose) cat(x)
  
  
  # Intersect sample names in study data and imprintome
  shared_sample_ids <- intersect(rownames(study_data), colnames(betas))
  study_data <- study_data[rownames(study_data) %in% shared_sample_ids, ]
  betas <- betas[, colnames(betas) %in% shared_sample_ids]
  verbosecat(sprintf("Shared sample IDs from betas and study data: %i \n", length(shared_sample_ids)))
  
  
  if (m_value_transform){
  t_meth = Matrix::t(sesame::BetaValueToMValue(betas)) %>% as.data.frame()
  } else {t_meth = Matrix::t(betas) %>% as.data.frame()  }
  
  # Parse model formula and assign the response and predictors
  out <- parse_lm_formula(model_str)
  
  if (out$response == "beta"){R = t_meth} else {R = study_data %>% dplyr::select(out$response)}
  if (out$primary_predictor == "beta"){P = t_meth} else {P =  study_data %>% dplyr::select(out$primary_predictor)}
  Pe =  study_data %>% dplyr::select(out$covariates)
  
  # If n_p_adj not specified, correct based on number of models.
  if(is.null(n_p_adj))  n_p_adj = max(c(ncol(R), ncol(P)))
  
  if (rm.na.all) rm.na.R <- rm.na.P <- rm.na.Pe <- rm.na.Pe <- TRUE
 
  
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
  
  # Remove all rows marked for removal
  rm.na.flags <- is.R.na | is.P.na | is.Pe.na
  verbosecat(sprintf(">> rm.na: Removing %.0f rows total (before imputation)...\n",
                     sum(rm.na.flags)))
  verbosecat(sprintf(">> rm.na:%.0f rows now remain.\n",
                     sum(!rm.na.flags)))
  R  <- R[!rm.na.flags, , drop = FALSE]
  P  <- P[!rm.na.flags, , drop = FALSE]
  Pe <- Pe[!rm.na.flags, , drop = FALSE]
  
  
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
                                             Pe = Pe,  family = family,
                                             impute_na = impute_na, .fit_model = fit_model,
                                             db_flag = db_flag, verbose = verbose)
  } else if (length(Pind) > 1) {
    # For each predictor variable
    foreach_fun <-  function (x) GLM_parallel(R = R, Rind = Rind, P = P, Pind = x,
                                              Pe = Pe,  family = family,
                                              impute_na = impute_na, .fit_model = fit_model,
                                              db_flag = db_flag, verbose = verbose)
  } else {
    foreach_fun <- function (x) GLM_parallel(R = R, Rind = Rind, P = P, Pind = Pind,
                                             Pe = Pe,  family = family,
                                             impute_na = impute_na, .fit_model = fit_model,
                                             db_flag = db_flag, verbose = verbose)
  }
  test_run <- foreach_fun(1)
  verbosecat(sprintf("Formula: %s\n", test_run$Formula))
  
  
  # Identify parallel index processing
  if (length(Rind) > 1) {
    par_ind <- Rind
  } else if (length(Pind) > 1) {
    par_ind <- Pind
  } else {
    par_ind <- 1L
  }
  
  # Determine number of workers
  if (is.null(n.cores)) {
    n.cores <- max(1L, parallel::detectCores() - 1L)
  }
  
  # Create cluster
  cl <- snow::makeCluster(n.cores)
  
  # Ensure that the cluster is stopped even if an error occurs
  on.exit(
    try(snow::stopCluster(cl), silent = TRUE),
    add = TRUE
  )
  
  doSNOW::registerDoSNOW(cl)
  
  # Initialize progress-bar objects
  pb <- NULL
  opts <- NULL
  
  if (verbose) {
    pb <- utils::txtProgressBar(
      min = 0,
      max = length(par_ind),
      style = 3
    )
    
    # Ensure that the progress bar is closed if an error occurs
    on.exit(
      try(close(pb), silent = TRUE),
      add = TRUE
    )
    
    opts <- list(
      progress = function(n) {
        utils::setTxtProgressBar(pb, n)
      }
    )
  }
  
  # Fit models
  df_fits <-
    foreach::foreach(
      x = par_ind,
      .combine = rbind,
      .export = c("fit_model", "GLM_parallel"),
      .packages = "magrittr",
      .options.snow = opts
    ) %dopar% {
      foreach_fun(x)
    }
  
  # Guarantee a data frame for downstream `$` operations
  df_fits <- as.data.frame(df_fits, stringsAsFactors = FALSE)
  
  if (verbose) {
    close(pb)
    pb <- NULL
    cat("\n")
  }
  
  snow::stopCluster(cl)
  cl <- NULL
  
  finish <- Sys.time()
  verbosecat(" Finished.\n")
  
  
  # Export for debugging
  if (db_flag) {save(list = ls(all.names = TRUE), file = "study_imprint2.RData")}
  # load(file = "study_imprint2.RData")
  
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
  # If cpg or ICR sites are a predictor, name the dataframe "imp_site"
  if (!is.null(P)) {
    df_temp <- list(df_fits[df_fits$Variable %in% colnames(P),])
    names(df_temp) <- "imp_site"
    dfs_sep <- c(df_temp, dfs_sep)
  }
  
  # Set NA p-values to a max value of 1 (for sorting)
  na_fun <- function(df) { df$P_VAL[is.na(df$P_VAL)] <- 1; return(df)}
  # Sort each variable by p-value
  sort_fun <- function(df) {df=df[order(df$P_VAL, decreasing = FALSE),]; return(df)}
  
  
  # Calculate adjusted p_value
  verbosecat("Adjusting p-values...\n")
  adj_p_val <- function(df) {
    df$ADJ_P_VAL <- apply_fdr_correction(pvals = df$P_VAL, n = n_p_adj)
    return(df)
  }
  dfs_sorted <- lapply(dfs_sep, na_fun)
  dfs_sorted <- lapply(dfs_sorted, sort_fun)
  dfs_corr   <- lapply(dfs_sorted, adj_p_val)
  dfs_corr$example_formula <- test_run$Formula[1]
  
  # Print out results of analysis
  try(expr = {
    summarize_study(dfs_corr, varnames = "imp_site", max_p_val, print_sites = FALSE,
                    print_confounders = FALSE)
  })
  
  
  return(dfs_corr)
}





#' analyze_association
#'
#' @description fits a series of general linear models with parallel processing
#' using dataframes of: R response variables, P predictor variables (split across 
#' models), and Pe predictor variables (consistent across models). The model 
#' formula takes on the form:
#'
#' If R has many columns and P has one, models are parallelized over R:
#'
#' R\[,i] ~ P\[,1] + Pe\[,1] + Pe\[,2] + Pe\[,3] ...
#'
#' If P has many variables/ columns and R has one, models are parallelized over P:
#'
#' R\[,1] ~ P\[,i] + Pe\[,1] + Pe\[,2] + Pe\[,3] ...
#'
#' This designed that either R or P are the CpG beta values or the ICR beta
#' values. The code handles either case. Note that R and P cannot both be
#' multicolumn.
#'
#' @param R dataframe of response variable so  patients (rows) x variables
#'  (columns) could either be beta values for cpg or ICR sites, or
#'  study metadata.
#' @param P parallelized predictor variable(s)- data frame that is patients
#'  (rows) x variables (columns), could either be beta values for cpg or ICR
#'  sites. For study metadata, use Pe input argument.
#' @param Pe extra predictor variables that are not parallelized (included in all
#'  fitted models).
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
#' @param rm.na.all TODO
#' @param print_confounders TODO
#' @param verbose TODO
#' @param icr_mapping TODO
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %dopar%
#' @importFrom rlang .data
#'
#' @return results with fitted coefficients from the glm, sorted by p-value
#' @export
analyze_association <- function (R, P = NULL, Pe = NULL, family, n_p_adj = max(c(ncol(R), ncol(P))),
                           max_p_val = 0.05, impute_na = TRUE, n.cores = NULL,
                           db_flag = FALSE, rm.na.R = FALSE, rm.na.P = FALSE,
                           rm.na.Pe = FALSE, rm.na.all = FALSE,
                           print_confounders = FALSE, verbose = TRUE,
                           icr_mapping = NULL) {
  if (db_flag) {save(list = ls(all.names = TRUE), file = "analyze_association.RData")}
  # load(file = "analyze_association.RData")

  
  warning("This functio is being deprecated, please switch to imprintome_glm.")
  
  if (rm.na.all) rm.na.R <- rm.na.P <- rm.na.Pe <- rm.na.Pe <- TRUE
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



  # Remove all rows marked for removal
  rm.na.flags <- is.R.na | is.P.na | is.Pe.na
  verbosecat(sprintf(">> rm.na: Removing %.0f rows total (before imputation)...\n",
              sum(rm.na.flags)))
  verbosecat(sprintf(">> rm.na:%.0f rows now remain.\n",
              sum(!rm.na.flags)))

  R  <- R[!rm.na.flags, , drop = FALSE]
  P  <- P[!rm.na.flags, , drop = FALSE]
  Pe <- Pe[!rm.na.flags, , drop = FALSE]



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
                                             Pe = Pe,  family = family,
                                             impute_na = impute_na, .fit_model = fit_model,
                                             db_flag = db_flag, verbose = verbose)
    # test_run <- GLM_parallel(R = R, Rind = 1, P = P, Pind = Pind,
    #                          Pe = Pe,  family = family,
    #                          impute_na = impute_na,
    #                          db_flag = db_flag, verbose = verbose)
    test_run <- foreach_fun(1)
  } else if (length(Pind) > 1) {
    # For each predictor variable
    foreach_fun <-  function (x) GLM_parallel(R = R, Rind = Rind, P = P, Pind = x,
                                              Pe = Pe,  family = family,
                                              impute_na = impute_na, .fit_model = fit_model,
                                              db_flag = db_flag, verbose = verbose)
    # test_run <- GLM_parallel(R = R, Rind = Rind, P = P, Pind = 1,
    #                          Pe = Pe,  family = family,
    #                          impute_na = impute_na,
    #                          db_flag = db_flag, verbose = verbose)
    test_run <- foreach_fun(1)
  } else {
    foreach_fun <- function (x) GLM_parallel(R = R, Rind = Rind, P = P, Pind = Pind,
                                             Pe = Pe,  family = family,
                                             impute_na = impute_na, .fit_model = fit_model,
                                             db_flag = db_flag, verbose = verbose)
    # test_run <- GLM_parallel(R = R, Rind = 1, P = P, Pind = 1,
    #                          Pe = Pe,  family = family,
    #                          impute_na = impute_na,
    #                          db_flag = db_flag, verbose = verbose)
    test_run <- foreach_fun(1)
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

  # https://support.bioconductor.org/p/132108/
  # BiocParallel::MulticoreParam(workers = n.cores)
  # bp_par <- BiocParallel::SnowParam(workers = n.cores, type = "SOCK", exportvariables = TRUE)
  # bp_out <- BiocParallel::bpoptions(exports = c("fit_model", "GLM_parallel"),packages = c("magrittr","stats"))
  # Add progress bar
  if (verbose) {
    pb <- utils::txtProgressBar(max = max(par_ind), style = 3)
    opts <- list(progress =  function(n) utils::setTxtProgressBar(pb, n))
  } else { opts = NULL }
    # s1 <- Sys.time()
    # par_out <- BiocParallel::bplapply(
    #   X = par_ind[1:100], GLM_parallel, R = R, Rind = Rind, P = P,
    #   Pe = Pe,  family = family,
    #   impute_na = impute_na, .fit_model = fit_model,
    #   db_flag = db_flag, verbose = FALSE, BPPARAM = bp_par, BPOPTIONS = bp_out)
    # # df_par <- do.call(rbind, par_out)
    # Sys.time()-s1 
    # 
  # Execute parallel processing                               ##################
  #_____________________________________________________________________________
  df_fits_list <- 
    foreach::foreach(
      x = par_ind, combine = rbind, .export = c("fit_model", "GLM_parallel"), 
      .packages = c("magrittr"), .options.snow = opts
    ) %dopar% {
      foreach_fun(x)
    }
  df_fits <- do.call(rbind, df_fits_list)
  # End processing time
  cat("\n");  finish <- Sys.time()
  verbosecat(" Finished.\n")


  # Export for debugging
  if (db_flag) {save(list = ls(all.names = TRUE), file = "study_imprint2.RData")}
  # load(file = "study_imprint2.RData")
  
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
  # If cpg or ICR sites are a predictor, name the dataframe "imp_site"
  if (!is.null(P)) {
    df_temp <- list(df_fits[df_fits$Variable %in% colnames(P),])
    names(df_temp) <- "imp_site"
    dfs_sep <- c(df_temp, dfs_sep)
  }

  # Set NA p-values to a max value of 1 (for sorting)
  na_fun <- function(df) { df$P_VAL[is.na(df$P_VAL)] <- 1; return(df)}
  # Sort each variable by p-value
  sort_fun <- function(df) {df=df[order(df$P_VAL, decreasing = FALSE),]; return(df)}


  # Calculate adjusted p_value
  verbosecat("Adjusting p-values...\n")
  adj_p_val <- function(df) {
    df$ADJ_P_VAL <- apply_fdr_correction(pvals = df$P_VAL, n = n_p_adj)
    return(df)
  }
  dfs_sorted <- lapply(dfs_sep, na_fun)
  dfs_sorted <- lapply(dfs_sorted, sort_fun)
  dfs_corr   <- lapply(dfs_sorted, adj_p_val)
  dfs_corr$example_formula <- test_run$Formula[1]

  # Print out results of analysis
  try(expr = {
    summarize_study(dfs_corr, varnames = "imp_site", max_p_val, print_sites = FALSE,
                    print_confounders = FALSE)
  })


  return(dfs_corr)
}



# TODO: fill in what this does
# lambda <-
#   stats::median(stats::qchisq( as.numeric( as.character(df_fits_sorted$P_VAL)),
#                                df = 1, lower.tail = F),
#                 na.rm = T) / stats::qchisq( 0.5, 1)
# verbosecat(sprintf("> Lambda: %f\n", lambda))



