


#' Fit Methylation Association Models Across CpG or ICR Sites
#'
#' Fits generalized linear models with methylation as either the response
#' or the first predictor, adjusting for additional study variables.
#'
#' @param model_str Character model formula using study_data column names
#'   and the special term beta for methylation. Put beta on the left-hand
#'   side or as the first term on the right-hand side, for example
#'   "outcome ~ beta + age" or "beta ~ exposure + age".
#' @param study_data Data frame with one row per sample. Row names must
#'   match sample identifiers in the columns of betas.
#' @param betas Numeric methylation matrix or data frame, sites in rows and
#'   samples in columns. Row names identify CpG or ICR sites.
#' @param family GLM family specification passed to the model fitter,
#'   usually "gaussian" or "binomial".
#' @param n_p_adj Number of comparisons for apply_fdr_correction(). NULL
#'   uses the larger number of response or primary-predictor columns.
#'   The adjustment helper requires this not to exceed the number of
#'   p-values in each adjusted result table.
#' @param max_p_val Adjusted p-value threshold for printed summaries.
#'   Does not filter the returned results.
#' @param impute_na Logical; enable MICE imputation in fit_model_glm() for
#'   missing model values remaining after the requested row removals.
#' @param n.cores Number of parallel workers. NULL uses the detected core
#'   count minus one, with a minimum of one.
#' @param db_flag Logical; save debugging environments to imprintome_glm.RData
#'   and study_imprint2.RData, and enable debugging in single_imprint_glm().
#' @param rm.na.R Logical; remove samples missing any response-column value
#'   before fitting any models.
#' @param rm.na.P Logical; remove samples missing any primary-predictor
#'   column value before fitting any models.
#' @param m_value_transform Logical; convert betas to M-values using
#'   sesame::BetaValueToMValue(). Defaults to TRUE.
#' @param rm.na.Pe Logical; remove samples missing any covariate value.
#' @param rm.na.all Logical; enable all three row-removal flags.
#' @param verbose Logical; print progress and model information. The final
#'   summary is attempted independently of this flag.
#'
#' @details
#' Only samples shared by study_data and betas are retained. Predictor rows
#' are aligned to response rows by sample ID. Models vary over methylation
#' sites while retaining the same study covariates. Requested missing-value
#' removal is performed across the entire selected response or predictor
#' table, so missingness at one site can exclude a sample from every model.
#' Remaining missingness is handled by fit_model_glm().
#'
#' Each result table is sorted by raw p-value, replacing missing raw p-values
#' with 1 before applying the package's FDR correction. Fitting failures caught
#' by single_imprint_glm() can produce placeholder rows with missing estimates and
#' p-values of 1; see that function for debug-file behavior.
#'
#' @return A named list of coefficient-result data frames. The imp_site
#'   element combines coefficients for primary predictors; other elements
#'   contain results for individual additional coefficient names. When beta
#'   is the response, methylation site IDs appear in Response; when beta is
#'   the first predictor, they appear in Variable.
#'   Tables contain Response, Variable, Estimate, StdError, Statistic, P_VAL,
#'   Confounder, ADJ_P_VAL, Family, Formula, Model_Id, and aic, as described
#'   in single_imprint_glm(), with ADJ_P_VAL filled by apply_fdr_correction().
#'   The example_formula attribute holds the first fitted formula.
#' @seealso [single_imprint_glm()], [fit_model_glm()], [apply_fdr_correction()]
#' @importFrom magrittr %>%
#' @importFrom foreach %dopar%
#' @importFrom rlang .data
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
  #_____________________________________________________________________________
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
  #-____________________________________________________________________________
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
  #_____________________________________________________________________________
  Rind <- 1: max(is.null(ncol(R)),ncol(R))
  Pind <- 1: max(is.null(ncol(P)),ncol(P))
  
  # Print short label for model (excluding confounders)
  verbosecat(sprintf("First model label: %s ~ %s ...\n", colnames(R)[1],
                     colnames(P)[1]))
  verbosecat("Example output of first model:\n")
  
  # Define parallel processing function for GLM
  #_____________________________________________________________________________
  # Run function on first input in verbose to show model for debugging
  if (length(Rind) > 1) {
    # For each response variable
    foreach_fun <- function (x) single_imprint_glm(R = R, Rind = x, P = P, Pind = Pind,
                                             Pe = Pe,  family = family,
                                             impute_na = impute_na, .fit_model = fit_model_glm,
                                             db_flag = db_flag, verbose = verbose)
  } else if (length(Pind) > 1) {
    # For each predictor variable
    foreach_fun <-  function (x) single_imprint_glm(R = R, Rind = Rind, P = P, Pind = x,
                                              Pe = Pe,  family = family,
                                              impute_na = impute_na, .fit_model = fit_model_glm,
                                              db_flag = db_flag, verbose = verbose)
  } else {
    foreach_fun <- function (x) single_imprint_glm(R = R, Rind = Rind, P = P, Pind = Pind,
                                             Pe = Pe,  family = family,
                                             impute_na = impute_na, .fit_model = fit_model_glm,
                                             db_flag = db_flag, verbose = verbose)
  }
  test_run <- foreach_fun(1)
  verbosecat(sprintf("Formula: %s\n", test_run$Formula[1]))
  
  
  # Identify parallel index processing
  if (length(Rind) > 1) {
    par_ind <- Rind
  } else if (length(Pind) > 1) {
    par_ind <- Pind
  } else {
    par_ind <- 1L
  }
  
  # Determine number of workers
  if (is.null(n.cores)) {n.cores <- max(1L, parallel::detectCores() - 1L)  }
  # Create cluster
  cl <- snow::makeCluster(n.cores)
  
  # Ensure that the cluster is stopped even if an error occurs
  on.exit(try(snow::stopCluster(cl), silent = TRUE),add = TRUE  )
  
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
    on.exit(try(close(pb), silent = TRUE), add = TRUE )
    
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
      .export = c("fit_model_glm", "single_imprint_glm"),
      .packages = "magrittr",
      .options.snow = opts
    ) %dopar% {
      foreach_fun(x)
    }
  
  # Guarantee a data frame for downstream `$` operations
  df_fits <- as.data.frame(df_fits, stringsAsFactors = FALSE)
  
  if (verbose) {close(pb); pb <- NULL; cat("\n")}
  
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
  
  # Go through each model variable and extract fit and stats
  dfs_sep <- list()
  for (n in seq_along(model_vars)) {
    dfs_sep[[n]] <- df_fits[df_fits$Variable == model_vars[n],]
  }
  names(dfs_sep) <- model_vars
  # If cpg or ICR sites are a predictor, label the dataframe "imp_site"
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
  attr(dfs_corr, "example_formula") <- test_run$Formula[1]

  
  # Print out results of analysis
  try(expr = {
    summarize_study(dfs_corr, varnames = "imp_site", max_p_val, print_sites = FALSE,
                    print_confounders = FALSE)
  })
  
  
  return(dfs_corr)
}






#' summarize_study
#'
#' @description Prints out results of a study analysis from analyze_association.
#' @param dfs a list of dataframes that give model output results for each
#' variable in model.
#' @param varnames vector of strings specifying which variables to cummarize.
#' @param max_p_val maximum p-value threshold for reported results. Both
#' adjusted and unadjusted p-values are reported.
#' @param print_sites boolean, when true, will print all cpg/icr sites that are
#' statistically significant (adjust p-value)
#' @param print_confounders boolean, when true, will print out statistical
#' summary of confoudner variables
#' @export
summarize_study <- function(dfs, varnames = NULL, max_p_val = 0.05,
                            print_sites = TRUE, print_confounders = FALSE) {
  cat(sprintf("Formula: %s \n", dfs$example_formula))
  
  if (is.null(varnames)) varnames <- names(dfs)
  
  sig_list = list()
  for (n in seq_along(varnames)) {
    if (dfs[[n]]$Confounder[1] == 0 || print_confounders ) {
      cat(sprintf("%s:\n", varnames[n]))
      cat(sprintf(">>  %.0f imprint sites have p_val < %.2f\n",
                  sum(dfs[[n]]$P_VAL < max_p_val), max_p_val))
      cat(sprintf(">>  %.0f imprint sites have adj_p_val < %.2f\n",
                  sum(dfs[[n]]$ADJ_P_VAL < max_p_val), max_p_val))
      if (sum(dfs[[n]]$ADJ_P_VAL < max_p_val)>0 && print_sites) {
        
        print(dfs[[n]][dfs[[n]]$ADJ_P_VAL < max_p_val,] %>%
                dplyr::select(,-c("Formula", "Model_Id")))
      }
      
      # Print ho wmany model fittings failed
      cat(sprintf("%.0f/ %.0f of model fits failed.\n", sum(is.na(dfs[[n]]$Estimate)), 
                  nrow(dfs[[n]])))
      cat("\n")
    }
    
    sig_list[[n]] <- dfs[[n]][ dfs[[n]]$ADJ_P_VAL < max_p_val, ]
  }
  
  df_sig = do.call(rbind, sig_list)
}





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
#' @param db_flag Logical; save the initial environment to single_imprint_glm.RData
#'   in the working directory. Failure snapshots can be written even if FALSE.
#' @param .fit_model Fitting function with arguments model_data,
#'   formula_string, family, impute_na, and n_imputes, in that order.
#'   Must return a list containing cf (a coefficient table including an
#'   intercept row) and aic. Defaults to fit_model_glm().
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
#' single_imprint_glm_error_pid-<process ID>.RData snapshot. A constant response
#' detected before fitting writes single_imprint_glm_error.RData and stops.
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
#' @seealso [fit_model_glm()], [imprintome_glm()]
#' @importFrom magrittr %>%
#' @export
single_imprint_glm = function(R, Rind = 1, P = NULL, Pind = 1, Pe = NULL,
                        family = "binomial", verbose = FALSE, impute_na = TRUE,
                        db_flag = FALSE, .fit_model = fit_model_glm) {
  if (db_flag) {save(list = ls(all.names = TRUE), file = "single_imprint_glm.RData") }
  # load(file = "single_imprint_glm.RData")
  
  stopifnot("single_imprint_glm:error: R must be a dataframe" = is.data.frame(R) ||
              is.matrix(R))
  stopifnot("single_imprint_glm:error: P must be a dataframe" = is.data.frame(P) ||
              is.matrix(P) || is.null(P))
  if (length(Rind)!=1 || length(Pind)!=1) {
    stop(sprintf(paste0("single_imprint_glm:error: both length() of Rind(==%.0f) and ",
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
    save(list = ls(all.names = TRUE), file = "single_imprint_glm_error.RData") # load(file = "single_imprint_glm.RData")
    stop(sprintf("single_imprint_glm: Response variable only has 1 unique value at
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
    save(list = ls(all.names = TRUE), file = paste0("single_imprint_glm_error_pid-", Sys.getpid(),".RData"))
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
#' @seealso [single_imprint_glm()]
#' @export
fit_model_glm <- function(model_data, formula_string, family, impute_na, n_imputes) {
  
  
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






#' Parse a Linear Model Formula
#'
#' Parses a model formula into its response variable, primary predictor, and
#' covariates. The first term on the right-hand side of the formula is treated
#' as the primary predictor, and all remaining right-hand-side terms are
#' returned as covariates.
#'
#' The input may be supplied either as a formula object or as a character
#' string that can be converted to a formula.
#'
#' @param formula A model formula or character string specifying a model.
#' The first term on the right-hand side is interpreted as the primary
#' predictor, while any subsequent terms are interpreted as covariates.
#'
#' @return A named list with the following elements:
#' \describe{
#' \item{response}{Character vector containing the variable name(s) appearing
#' in the response expression.}
#' \item{primary_predictor}{Character string containing the first
#' right-hand-side term, or \code{NULL} if no right-hand-side terms are
#' present.}
#' \item{covariates}{Character vector containing all right-hand-side terms
#' after the primary predictor. Returns \code{character(0)} if no
#' covariates are present.}
#' }
#'
#' @details
#' Right-hand-side terms are obtained from \code{stats::terms()}. Consequently,
#' compound terms such as interactions are retained using their model-term
#' representation (for example, \code{"x"}).
#'
#' The response is extracted using \code{all.vars()}, so a transformed response
#' such as \code{log(y)} is returned as \code{"y"} rather than
#' \code{"log(y)"}.
#'
#' @examples
#' parse_lm_formula(y ~ x + age + sex)
#'
#' parse_lm_formula("y ~ treatment + age + sex")
#'
#' parse_lm_formula(y ~ x)
#'
#' @export
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
    covariates = covariates,
    all = c(response, primary_predictor, covariates)
  )
}
