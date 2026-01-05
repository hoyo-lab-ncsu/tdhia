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
#' (for a binomial outcome) or an F test (for a continuous outcome).

#' @param cpg_beta cpg beta matrix (cpg_id as COLUMNS x sample_id as ROWS).
#'  Values are assumed to be on beta scale (0-1).
#'  Use m_value_transform to convert to m-values (recommended).

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

#' @param icr_ids vector of strings of icr_ids to be tested. Default = NULL,
#' tests all icrs that are covered by the input cpg beta matrix.

#' @param min_cpg minimum number of cpg sites for an icr to be included in
#' results (default is 3).

#' @returns test
#' @export
#' @author Kate Everly

#' Note that there is a line of code to remove any patients that are missing
#' data (outcome or covariates). I recommend manually filtering that out first
#' so that you know your exact sample size prior to running (just my personal
#' preference)
#' Another quick check prior to starting is to ensure both the cpg_beta and
#' df_study have the patients ids in rows
library(clusterSim)
library(dplyr)
library(tibble)
library(tidyr)



pc_regression_test <- function (cpg_beta, data_norm_type="n1", pct_variance = 0.80, df_study,
                                outcome, covariates, Patient_ID, family,
                                icr_ids = NULL, min_cpg = 3, verbose = T, n.cores){
  verbosecat <- function(x) if (verbose) cat(x)
  icr_mapping = tdhia::mapping_cpg_icr_ids # load in map!
  cpg_mapping <- icr_mapping %>% # will select just the CpG and ICR ids from the icr_mapping file
    dplyr::select(CpG_id, ICR_id)

  # if icr_ids no supplied, scan for all icr_ids covered with cpg_ids
  if (is.null(icr_ids)) {
    icr_ids <- cpg_mapping %>%
      dplyr::filter(CpG_id %in% colnames(cpg_beta)) %>% dplyr::pull(ICR_id) %>%
      unname() %>% unique()
  }

  # Tranform to m-values if requested
  # need to work on this.... not sure if the CpG beta needs to get t first
  #if (m_value_transform) tZ = sesame::BetaValueToMValue(tZ)

  verbosecat("> Remove rows in df_study that contain NA...\n")
  # Remove all samples with NA values from study data
  df_study <- df_study %>% select(all_of(c(outcome, covariates, Patient_ID))) %>% tidyr::drop_na()

  # Make the samples match and order them the same
  verbosecat("> Forcing sample id order to match between data and df_study.\n")
  shared_sample_ids <- intersect(rownames(cpg_beta),  rownames(df_study))
  cpg_beta <- cpg_beta[shared_sample_ids, ]
  df_study <- df_study[shared_sample_ids, , drop = FALSE]
  #df_study <- df_study[rownames(df_study) %in% shared_sample_ids, ]
  # Order df_study the same sample order as cpg_beta
  df_study <- df_study[match(x = rownames(cpg_beta), table = rownames(df_study)),]

  # Dataframe to store results of test
  df_results <- data.frame(ICR_id = icr_ids, p_value = NA, n_cpg = NA)

  # packaged into function for single ICR to make parallel conversion easier in future
  if (n.cores==1){

    verbosecat("> PC regression processing on single core.\n")
    out <-  list()

    for (n in 1:length(icr_ids)) {

      # Get the CpG IDs of the CpGs for an ICR
      subset_cpg_ids <- cpg_mapping %>%
        dplyr::filter(ICR_id == icr_ids[n]) %>% #  dplyr::filter(.data$icr_id == .env$icr_id)
        pull(CpG_id)

      # Subset the beta matrix to only have data from the CpGs found in the ICR
      subset_cpg_beta <- cpg_beta %>%
        dplyr::select(any_of(subset_cpg_ids))

      # Data normalization prior to creating principal components
      subset_cpg_beta.norm <- data.Normalization(subset_cpg_beta , type=data_norm_type, normalization="column")

      # Create principal components
      icr.pca <- prcomp(subset_cpg_beta.norm, center=TRUE, scale.= TRUE)
      eigenvalues <- icr.pca$sdev^2 #Calculate eigenvalues
      proportion_variance <- eigenvalues / sum(eigenvalues) #calculate the proportion variance of each eigenvalue
      cumulative_variance <- cumsum(proportion_variance) #calculate the cumulative variance
      # Create a dataframe describing each PC
      pc_info <- tibble(
        PC = paste0("PC", seq_along(eigenvalues)),
        eigenvalue = eigenvalues,
        prop_var = proportion_variance,
        cum_var = cumulative_variance)
      # Select only the top PCs contributing to ____% variance
      if (nrow(pc_info) > 1) {
        # Find the last PC where cum_var is at least ____%
        pcs_pct_variance <- which(pc_info$cum_var >= pct_variance)
        if (length(pcs_pct_variance) > 0) {
          cutoff_index <- min(pcs_pct_variance)
          selected_pcs <- pc_info[1:cutoff_index, ]
        }
      } else if (nrow(pc_info) == 1) {
        selected_pcs <- pc_info
      }

      # Get the actual PC scores (samples × PCs), then keep only top PCs
      pcs <- as.data.frame(icr.pca$x) %>%
        dplyr::select(all_of(selected_pcs$PC))
      predictor_full <- pcs

# FOR BINOMIAL OUTCOMES
      if (family == "binomial") {

        predictor_cpg_names <- colnames(predictor_full) # define before adding the Patient_ID column or else Patient_ID will be considered a variable
        predictor_full <- predictor_full %>%
          mutate(Patient_ID = rownames(predictor_full))

        study_data_full <- df_study
        study_data_full[[outcome]] <- as.factor(study_data_full[[outcome]])
        print(names(study_data_full))
        #the below line will put ALL the predictors together (PCs and covariates)
        combined_data_full <- merge(study_data_full, predictor_full, by = "Patient_ID")
        rownames(combined_data_full) <- combined_data_full$Patient_ID
        combined_data_full$Patient_ID <- NULL

        formula_str_full <- paste(outcome, "~", paste(c(predictor_cpg_names, covariates), collapse = " + "))
        model_formula_full <- as.formula(formula_str_full)

        model_full <- tryCatch({
          withCallingHandlers(
            expr = {
              glm(model_formula_full, data = combined_data_full, family = "binomial")
            },
            warning = function(w) {
              message(sprintf("Warning in full model for ICR %s: %s", icr_ids[n], conditionMessage(w)))
              invokeRestart("muffleWarning")
            }
          )
        }, error = function(e) {
          message(sprintf("Error in full model for ICR %s: %s", icr_ids[n], conditionMessage(e)))
          return(NULL)  # No invokeRestart for errors
        })

        formula_str_red <- paste(outcome, "~", paste(covariates, collapse = " + "))
        model_formula_red <- as.formula(formula_str_red)

        model_red <- tryCatch({
          glm(model_formula_red, data = combined_data_full, family = "binomial")
        }, warning = function(w) {
          message(sprintf("Warning in reduced model for ICR %s: %s", icr_ids[n], conditionMessage(w)))
          invokeRestart("muffleWarning")
        }, error = function(e) {
          message(sprintf("Error in reduced model for ICR %s: %s", icr_ids[n], conditionMessage(e)))
          return(NULL)
        })

        #calculate likelihood ratio test statistic with anova and ChiSq test
        lrt <- anova(model_red, model_full, test = "Chisq")
        lrt_row <- tail(lrt, 1) #the last row has the test stats

        out[[n]] <- data.frame(
          ICR_id = icr_ids[n],
          lrt_p_value = lrt_row$`Pr(>Chi)`,
          n_cpg = length(subset_cpg_ids)
        )
      }


# FOR CONTINUOUS OUTCOMES
      if (family == "continuous") {
        predictor_cpg_names <- colnames(predictor_full) # define before adding the Patient_ID column or else Patient_ID will be considered a variable
        predictor_full <- predictor_full %>%
          mutate(Patient_ID = rownames(predictor_full))

        study_data_full <- df_study
        #study_data_full[[outcome]] <- as.factor(study_data_full[[outcome]])
        #the below line will put ALL the predictors together (PCs and covariates)
        combined_data_full <- merge(study_data_full, predictor_full, by = "Patient_ID")
        rownames(combined_data_full) <- combined_data_full$Patient_ID
        combined_data_full$Patient_ID <- NULL

        formula_str_full <- paste(outcome, "~", paste(c(predictor_cpg_names, covariates), collapse = " + "))
        model_formula_full <- as.formula(formula_str_full)

        model_full <- tryCatch({
          withCallingHandlers(
            expr = {
              lm(model_formula_full, data = combined_data_full)
            },
            warning = function(w) {
              message(sprintf("Warning in full model for ICR %s: %s", icr_ids[n], conditionMessage(w)))
              invokeRestart("muffleWarning")
            }
          )
        }, error = function(e) {
          message(sprintf("Error in full model for ICR %s: %s", icr_ids[n], conditionMessage(e)))
          return(NULL)  # No invokeRestart for errors
        })

        formula_str_red <- paste(outcome, "~", paste(covariates, collapse = " + "))
        model_formula_red <- as.formula(formula_str_red)

        model_red <- tryCatch({
          lm(model_formula_red, data = combined_data_full)
        }, warning = function(w) {
          message(sprintf("Warning in reduced model for ICR %s: %s", icr_ids[n], conditionMessage(w)))
          invokeRestart("muffleWarning")
        }, error = function(e) {
          message(sprintf("Error in reduced model for ICR %s: %s", icr_ids[n], conditionMessage(e)))
          return(NULL)
        })

        #calculate F test with ANOVA
        f_test <- anova(model_red, model_full)
        f_test_row <- tail(f_test, 1) #the last row has the test stats

        out[[n]] <- data.frame(
          ICR_id = icr_ids[n],
          p_value = f_test_row$`Pr(>F)`,
          n_cpg = length(subset_cpg_ids)
        )
      }


    }

# REGARDLESS OF BINOMIAL OR CONTINUOUS OUTCOME, COMBINE RESULTS AND CALCUALTE ADJ-P/Q VALUES
    df_results = do.call(rbind, out)

  } else {
    verbosecat("Ahhhhhh more than one core!\n
               Go back and change to one core until I've got multicore setup\n")
  }

  # Calculate adjusted p-value and q-value
  df_results$adj_p_value <- p.adjust(p = df_results$p_value, method = "fdr")
  df_results$q_value <- qvalue::qvalue(p = df_results$p_value, fdr.level = 0.05)$qvalues
  df_results <- df_results %>% dplyr::arrange(adj_p_value)

  verbosecat(sprintf("> Filtering out ICRs with < %d cpg sites...\n", min_cpg))
  df_results <- df_results%>% filter(n_cpg >= min_cpg)

  return(df_results)
}

