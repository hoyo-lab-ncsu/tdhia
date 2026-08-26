#' pc_regression_test
#'
#' Principal component regression followed by likelihood ratio test or F test
#' to test for the association of ICR methylation with an outcome. Both
#' continuous and dichotomous response variables are supported.
#' For each ICR, the CpG sites will be reduced to principal components that
#' account for a chosen amount of variance (default 80%). Two regression models
#' will be fit: a reduced model with only covariates and a full model with
#' covariates and principal components (ICR methylation). In the case of
#' binomial outcomes, logistic regression models using the binomial family are
#' fit. In the case of continuous outcomes, linear regression models are fit.
#' The full and reduced models are then compared with a likelihood ratio test
#' (for a binomial outcome) or an F test (for a continuous outcome). The resulting
#' q-value adjusts for multiple test correction and tells you if the full
#' model was significantly improved over the reduced.
#'
#' Note that there is a line of code to remove any patients that are missing
#' data (outcome or covariates). I recommend manually filtering that out first
#' so that you know your exact sample size prior to running (just my personal
#' preference). Another quick check prior to starting is to ensure both the cpg_beta and
#' df_study have the patients ids in rows.

#' @param cpg_beta cpg beta matrix (cpg_id as COLUMNS x sample_id as ROWS).
#'  Values are assumed to be on beta scale (0-1).
#'  Use m_value_transform to convert to m-values (recommended).

#' @param m_value_transform boolean, when true, transform beta values to
#' m-values to control heteroskedasticity. Default = TRUE.

#' @param data_norm_type string. Type of data normalization you want to perform
#' on the cpg beta values prior to creating the principal components. Options
#' are in the data.Normalization function of the clusterSim package
#'    Default is n1 (standardization)

#' @param pct_variance string, the minimum chosen level of variance that the
#' principal components will account for. Default is 80%.
#'  For example, if PC1 accounts for 59% of the cumulative variance, PC2 for
#'  72% of the of the cumulative variance, PC3 for 81% of the of the cumulative
#'  variance, and PC4 for 88% of the cumulative variance, the PCs that will
#'  be used for analysis are PC1, PC2 and PC3.
#'
#' @param df_study dataframe of sample associated data to be used in linear
#' models, (nrows = sample size). Columns should include those specified with
#' 'outcome' and 'covariates' input arguments. There should be a column labeled
#' "Patient_ID" to link the cpg_beta and df_study

#' @param outcome string, column name of response variable located in df_study.

#' @param covariates vector of strings of column names of predictors that are
#' located in df_study.
#'
#' @param Patient_ID string, whatever your patient ID column is called in your
#' study data

#' @param family string. Define your outcome type (categorical outcomes- ex.
#' case vs control) are defined as binomial. Continuous numerical outcomes as
#' continuous. This will determine if logistic regression with likelihood
#' ratio test, or linear regression with F test will be performed.
#'    "binomial" (default)
#'    "continuous"

#' @param icr_ids vector of strings of icr_ids to be tested. Default = NULL
#' (tests all icrs that are covered by the input cpg beta matrix).

#' @param min_cpg minimum number of cpg sites for an icr to be included in
#' results (default is 3).

#' @param n.cores number of cores to run on (default is 1)
#' @param verbose boolean, when true, prints additional output.
#' @returns test
#' @export
#' @author Kate Everly


pc_regression_test <- function (cpg_beta,
                                m_value_transform= TRUE,
                                data_norm_type="n1",
                                pct_variance = 0.80,
                                df_study,
                                outcome,
                                covariates,
                                Patient_ID,
                                family = "binomial",
                                icr_ids = NULL,
                                min_cpg = 3,
                                verbose = TRUE,
                                n.cores = 1){
  verbosecat <- function(x) if (verbose) cat(x)
  icr_mapping = tdhia::mapping_cpg_icr_ids # load in map!
  cpg_mapping <- icr_mapping %>%
    dplyr::select(CpG_id, ICR_id) %>% # will select just the CpG and ICR ids from the icr_mapping file
    dplyr::distinct()

  # if icr_ids no supplied, scan for all icr_ids covered with cpg_ids
  if (is.null(icr_ids)) {
    icr_ids <- cpg_mapping %>%
      dplyr::filter(CpG_id %in% colnames(cpg_beta)) %>% dplyr::pull(ICR_id) %>%
      unname() %>% unique()
  }

  # Tranform cpg_beta matrix from beta values to m-values if requested
  # even if you transform to m-values, the dataframe of cpg methylation levels is
  # still referred to as the "cpg_beta" throughout the script
  if (m_value_transform){
    cpg_beta <- sesame::BetaValueToMValue(cpg_beta)
  }

  verbosecat("> Remove rows in df_study that contain NA...\n")
  # Remove all samples with NA values from study data
  df_study <- df_study %>% dplyr::select(dplyr::all_of(c(outcome, covariates, Patient_ID))) %>% tidyr::drop_na()

  # Make the samples match and order them the same
  verbosecat("> Forcing sample id order to match between cpg_beta and df_study.\n")
  shared_sample_ids <- intersect(rownames(cpg_beta),  rownames(df_study))
  cpg_beta <- cpg_beta[shared_sample_ids, ]
  df_study <- df_study[shared_sample_ids, , drop = FALSE]
  df_study <- df_study[match(x = rownames(cpg_beta), table = rownames(df_study)),]

  # Safe_fit function is to catch errors with model fitting
  safe_fit <- function(expr, icr_id, model_type) {
    tryCatch(
      withCallingHandlers(
        expr = expr,
        warning = function(w) {
          message(sprintf("Warning in %s model for ICR %s: %s",
                          model_type, icr_id, conditionMessage(w)))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        message(sprintf("Error in %s model for ICR %s: %s",
                        model_type, icr_id, conditionMessage(e)))
        return(NULL)
      }
    )
  }

  # -------------------------------------------------------------------------
  # Helper function that processes ONE ICR (PCA + regression)
  # -------------------------------------------------------------------------

  process_icr <- function(icr_id){
    # Get the CpG IDs of the CpGs for an ICR
    subset_cpg_ids <- cpg_mapping %>%
      dplyr::filter(ICR_id == icr_id, CpG_id %in% colnames(cpg_beta)) %>%
      dplyr::pull(CpG_id) %>%
      unique()

    # Subset the beta matrix to only have data from the CpGs found in the ICR
    subset_cpg_beta <- cpg_beta %>%
      dplyr::select(dplyr::any_of(subset_cpg_ids))

    # Data normalization prior to creating principal components
    subset_cpg_beta_norm <- clusterSim::data.Normalization(subset_cpg_beta , type=data_norm_type, normalization="column")

    # Create Principal Components
    # -------------------------------------------------------------------------
    icr.pca <- stats::prcomp(subset_cpg_beta_norm, center = TRUE, scale. = TRUE)
    eigenvalues <- icr.pca$sdev^2
    prop_var <- eigenvalues / sum(eigenvalues) #calculate the proportion variance of each eigenvalue
    cum_var <- cumsum(prop_var) #calculate the cumulative variance
    # Select PCs contributing to pct_variance
    if (length(cum_var) > 1) {
      cutoff_index <- which(cum_var >= pct_variance)[1]
    } else {
      cutoff_index <- 1
    }
    selected_pcs <- paste0("PC", seq_len(cutoff_index))
    # Get PC scores
    pcs <- as.data.frame(icr.pca$x) %>%
      dplyr::select(dplyr::all_of(selected_pcs))
    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------

    # Merge PCs with study data
    #rownames(pcs) <- pcs[[Patient_ID]]
    #pcs[[Patient_ID]] <- NULL
    #combined_data <- cbind(df_study, pcs[rownames(df_study), , drop = FALSE])
    # HARD ENFORCE alignment

    pcs[[Patient_ID]] <- rownames(pcs)
    combined_data <- dplyr::inner_join(
      df_study,
      pcs,
      by = Patient_ID
    )

    stopifnot(
      nrow(combined_data) == nrow(df_study),
      nrow(combined_data) == nrow(pcs),
      !anyDuplicated(df_study[[Patient_ID]]),
      !anyDuplicated(pcs[[Patient_ID]])
    )
    rownames(combined_data) <- combined_data[[Patient_ID]]
    combined_data[[Patient_ID]] <- NULL

    predictor_cpg_names <- selected_pcs

    # Build formulas
    # -------------------------------------------------------------------------
    full_formula <- stats::as.formula(
      paste(outcome, "~", paste(c(predictor_cpg_names, covariates), collapse = " + "))
    )
    red_formula <- stats::as.formula(
      paste(outcome, "~", paste(covariates, collapse = " + "))
    )
    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------


    # Fit models
    # -------------------------------------------------------------------------
    if (family == "binomial") {
      combined_data[[outcome]] <- as.factor(combined_data[[outcome]])
      model_full <- safe_fit(stats::glm(full_formula, data = combined_data, family = "binomial"),
                             icr_id, "full")
      model_red  <- safe_fit(stats::glm(red_formula, data = combined_data, family = "binomial"),
                             icr_id, "reduced")

      lrt <- stats::anova(model_red, model_full, test = "Chisq")
      pval <- utils::tail(lrt, 1)$`Pr(>Chi)`
    } else {
      model_full <- safe_fit(stats::lm(full_formula, data = combined_data),
                             icr_id, "full")
      model_red  <- safe_fit(stats::lm(red_formula,  data = combined_data),
                             icr_id, "reduced")
      ftest <- stats::anova(model_red, model_full)
      pval <- utils::tail(ftest, 1)$`Pr(>F)`
    }
    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------

    # SAVE RESULTS
    tibble::tibble(
      ICR_id = icr_id,
      p_value = pval,
      n_cpg = ncol(subset_cpg_beta)
    )
  }


  # -------------------------------------------------------------------------
  # Chose to execute on single core or multicore
  # -------------------------------------------------------------------------
  if (n.cores == 1) {

    verbosecat("> PC regression processing on single core.\n")
    out <- lapply(icr_ids, process_icr)

  } else {

    verbosecat("> PC regression processing on multiple cores.\n")

    # macOS/Linux → MulticoreParam
    # Windows → SnowParam
    if (Sys.info()[["sysname"]] == "Windows") {
      param <- BiocParallel::SnowParam(workers = n.cores)
    } else {
      param <- BiocParallel::MulticoreParam(workers = n.cores)
    }

    out <- BiocParallel::bplapply(icr_ids, process_icr, BPPARAM = param)
  }


  # REGARDLESS OF BINOMIAL OR CONTINUOUS OUTCOME, COMBINE RESULTS AND CALCUALTE ADJ-P/Q VALUES
  df_results = do.call(rbind, out)
  df_results$adj_p_value <- p.adjust(p = df_results$p_value, method = "fdr")
  df_results$q_value <- qvalue::qvalue(p = df_results$p_value)$qvalues
  df_results <- df_results %>% dplyr::arrange(adj_p_value)

  verbosecat(sprintf("> Filtering out ICRs with < %d cpg sites...\n", min_cpg))
  df_results <- df_results%>% dplyr::filter(n_cpg >= min_cpg)

  return(df_results)
}
