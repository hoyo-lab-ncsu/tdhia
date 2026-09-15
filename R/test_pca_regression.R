#' Test ICR Associations Using Principal Component Regression
#'
#' Reduces each ICR's CpG measurements to principal components, then compares
#' a covariate-only model with a model containing covariates and those
#' components. Uses a likelihood-ratio test for binomial outcomes and an
#' F test for continuous outcomes.
#'
#' @param cpg_beta Numeric data frame of beta values with samples in rows and
#'   CpG IDs in columns. Row names must match sample IDs in df_study.
#' @param m_value_transform Logical; convert beta values to M-values before
#'   normalization and PCA. Defaults to TRUE.
#' @param data_norm_type Normalization type passed to
#'   clusterSim::data.Normalization(); "n1" standardizes columns.
#' @param pct_variance Numeric cumulative variance fraction, normally in
#'   (0, 1], used to select components separately for each ICR. Default: 0.80.
#' @param n_pcs Optional positive integer overriding pct_variance. Uses the
#'   requested number of components, capped at the number available per ICR.
#' @param df_study Data frame containing outcome, covariates, and sample IDs.
#'   Row names must match cpg_beta; ID-column values must agree with row names.
#' @param outcome Character string naming the response column in df_study.
#' @param covariates Character vector of adjustment-column names in df_study.
#'   The current reduced-model formula construction requires covariates.
#' @param Patient_ID Character string naming the sample-ID column selected
#'   from df_study. Use "Patient_ID": the per-ICR helper currently joins on
#'   that literal column name regardless of this argument.
#' @param family Character string. "binomial" fits logistic models and uses a
#'   chi-squared likelihood-ratio test. Other values, including "gaussian" and
#'   "continuous", use linear models and an F test. There is no default.
#' @param icr_ids Character vector of ICR IDs, or NULL to use ICRs represented
#'   by the CpG columns in the package mapping.
#' @param min_cpg Minimum CpG count for retaining an ICR in the final table.
#'   This filter is applied after fitting and multiple-testing adjustment.
#' @param verbose Logical; print progress. Model warnings and errors can
#'   still be reported when FALSE.
#' @param n.cores Number of workers. One uses lapply(); larger values use
#'   BiocParallel with SnowParam on Windows and MulticoreParam elsewhere.
#' @param db_flag Logical; save the initial environment to
#'   pc_regression_test.RData in the working directory. Defaults to FALSE.
#'
#' @details
#' Study rows with missing outcome, covariates, or sample IDs are removed,
#' then study and methylation data are aligned by row names. CpG measurements
#' are not imputed here. Each ICR is normalized by column, then prcomp() is
#' called with centering and scaling. Fitting warnings are reported; caught
#' fitting errors return NULL to subsequent model-comparison code and do not
#' guarantee that processing of other ICRs will continue.
#'
#' @return A table sorted by adj_p_value, with:
#'   - ICR_id: region identifier.
#'   - raw_p_value: full-versus-reduced model comparison p-value.
#'   - n_cpg and n_pc: numbers of CpGs and retained components.
#'   - model_converged: the full model's convergence field when supplied by
#'     the model object; lm objects do not provide this field.
#'   - model_range: range of fitted values from the full model.
#'   - adj_p_value: Benjamini-Hochberg adjusted p-value.
#'   - q_value: q-value calculated by qvalue::qvalue().
#'   Adjustments include fitted ICRs subsequently removed by min_cpg.
#' @seealso [skat_icr_test()], [tdhia_stat_tests()]
#' @author Kate Everly
#' @export
pc_regression_test <- function (
    cpg_beta, m_value_transform = TRUE,  data_norm_type="n1",
    pct_variance = 0.80, n_pcs = NULL,
    df_study, outcome, covariates, Patient_ID, family, icr_ids = NULL,
    min_cpg = 3, verbose = TRUE, n.cores = 1, db_flag = FALSE) {
  
  if (db_flag) {save(list = ls(all.names = TRUE), file = "pc_regression_test.RData")}
  # load(file = "pc_regression_test.RData")
  
  verbosecat <- function(x) if (verbose) cat(x)
  icr_mapping = tdhia::mapping_cpg_icr_ids # load in map!
  cpg_mapping <- icr_mapping %>%
    dplyr::select(CpG_id, ICR_id) %>% # will select just the CpG and ICR ids from the icr_mapping file
    dplyr::distinct()

  
  if (!is.null(n_pcs)) {
    if (length(n_pcs) != 1L || !is.numeric(n_pcs) || is.na(n_pcs) ||
        !is.finite(n_pcs) || n_pcs < 1 ||  n_pcs != as.integer(n_pcs)) {
      stop("`n_pcs` must be NULL or a single positive integer.")
    }
    n_pcs <- as.integer(n_pcs)
  }
  
  
  # if icr_ids no supplied, scan for all icr_ids covered with cpg_ids
  if (is.null(icr_ids)) {
    icr_ids <- cpg_mapping %>%
      dplyr::filter(CpG_id %in% colnames(cpg_beta)) %>% dplyr::pull(ICR_id) %>%
      unname() %>% unique()
  }

  # Transform cpg_beta matrix from beta values to m-values if requested
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

  # -------------------------------------------------------------------------
  # Chose to execute on single core or multicore
  # -------------------------------------------------------------------------
  if (n.cores == 1) {
    verbosecat("> PC regression processing on single core.\n")
    out <- lapply(
      X = icr_ids, FUN =  pcr_single_icr, cpg_beta = cpg_beta, df_icr_pcr = df_icr_pcr,
      data_norm_type = data_norm_type, cpg_mapping = cpg_mapping, n_pcs = n_pcs, 
      df_study = df_study, outcome = outcome, covariates= covariates, safe_fit = safe_fit,
      pct_variance = pct_variance, family = family)

  } else {
    verbosecat("> PC regression processing on multiple cores.\n")
    # macOS/Linux → MulticoreParam
    # Windows → SnowParam
    if (Sys.info()[["sysname"]] == "Windows") {
      param <- BiocParallel::SnowParam(workers = n.cores)
    } else {
      param <- BiocParallel::MulticoreParam(workers = n.cores)
    }
    out <- BiocParallel::bplapply(
      X = icr_ids, FUN = pcr_single_icr, cpg_beta = cpg_beta,df_icr_pcr = df_icr_pcr,
      data_norm_type = data_norm_type, cpg_mapping = cpg_mapping, n_pcs = n_pcs, 
      df_study= df_study,outcome = outcome,covariates = covariates, safe_fit = safe_fit,
      pct_variance = pct_variance, family = family, BPPARAM = param)
  }


  # REGARDLESS OF BINOMIAL OR CONTINUOUS OUTCOME, COMBINE RESULTS AND CALCUALTE ADJ-P/Q VALUES
  df_results = do.call(rbind, out)
  df_results$adj_p_value <- p.adjust(p = df_results$raw_p_value, method = "fdr")
  df_results$q_value <- qvalue::qvalue(p = df_results$raw_p_value, fdr.level = 0.05)$qvalues
  df_results <- df_results %>% dplyr::arrange(adj_p_value)

  verbosecat(sprintf("> Filtering out ICRs with < %d cpg sites...\n", min_cpg))
  df_results <- df_results%>% dplyr::filter(n_cpg >= min_cpg)

  return(df_results)
}




pcr_single_icr <- function(icr_id, cpg_beta, df_icr_pcr, data_norm_type, cpg_mapping, n_pcs,df_study, outcome, covariates,safe_fit,pct_variance, family) {
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
  
  # Select the number of PCs based on pct_variance or n_pcs
  n_available_pcs <- ncol(icr.pca$x)
  if (is.null(n_pcs)) {
    # Default behavior: retain enough PCs to explain pct_variance
    cutoff_index <- which(cum_var >= pct_variance)[1]
  } else {
    # Override behavior: retain exactly n_pcs, subject to availability
    cutoff_index <- min(n_pcs, n_available_pcs)
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
  
  pcs$Patient_ID <- rownames(pcs)
  combined_data <- dplyr::inner_join(
    df_study,
    pcs,
    by = "Patient_ID"
  )
  
  stopifnot(nrow(combined_data) == nrow(df_study) || nrow(combined_data) == nrow(pcs))
  rownames(combined_data) <- combined_data$Patient_ID
  combined_data$Patient_ID <- NULL
  
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
  tibble::tibble(ICR_id = icr_id, raw_p_value = pval, n_cpg = ncol(subset_cpg_beta),
                 n_pc = cutoff_index, model_converged = model_full$converged, model_range = diff(range(fitted(model_full))))
}

