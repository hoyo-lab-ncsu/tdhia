


#' tdha_stat_tests
#' 
#' 
#' @param model_str todo
#' @param beta_data todo
#' @param study_data todo
#' @param transform_to_m todo
#' @param model_prefix todo
#' @param impute_na todo
#' @param n.cores todo
#' 
#' @export
tdhia_stat_tests <- function(
    model_str, betas, study_data, m_value_transform = T, data_cache_path, 
    model_prefix = "", impute_na = TRUE, n.cores = max(c(parallel::detectCores()-4, 1)),
    verbose = T, family = NULL,db_flag= T, overwrite_cache= F) {
  
  if (db_flag) {save(list = ls(all.names = TRUE), file = "tdhia_stat_tests.RData")}
  # load(file = "tdhia_stat_tests.RData")
  overwrite_cache= F
  
  verbosecat <- \(x) if (verbose) cat(x)
  
  # Parse model variables
  parsed_model <- tdhia::parse_lm_formula(model_str)
  

  # Extract patient id in rownames
  # study_data <- study_data %>% rownames_to_column("patient_id")
  
  # Extract columns used in study data
  study_columns <- c(parsed_model$all[parsed_model$all != "beta"])#,"patient_id")
  
  
  # Create beta data and study_data with NAs                                 ####
  #_____________________________________________________________________________
  # Set p_val threshold to NA
  study_data_na <- study_data %>% select(all_of(study_columns))
  betas_na <- betas
  betas_na$input_args$probe_data_cache <- betas_na$probe_beta
  betas_na$input_args$set_failed_betas_na <- T
  
  betas_na$input_args$db_flag=T
  betas_na <- do.call(what = tdhia::tdhia_pipeline, args = betas_na$input_args)
  # Discard patients that were filtered in processing
  study_data_na <-  study_data_na %>% rownames_to_column("patient_id") %>% 
    filter(patient_id %in% colnames(betas_na$cpg_beta$cpg_beta_df))
  rownames(study_data_na) <- study_data_na$patient_id
  
  # Create complete data with no NAs                                        ####
  #_____________________________________________________________________________
  study_data_complete <- study_data %>% select(all_of(study_columns)) %>% na.omit()
  betas_complete <- betas
  # Only include patient columns that are found in complete records in study data
  # betas_complete$input_args$probe_data_cache <- subset_probe_beta(data_beta$probe_beta, rownames())
  # Only load a subset
  betas_complete$input_args$idat_basenames <- rownames(study_data_complete)
  # Keeps all measurements
  betas_complete$input_args$set_failed_betas_na <- F
  betas_complete$input_args$db_flag=T
  betas_complete <- do.call(what = tdhia::tdhia_pipeline, args = betas_complete$input_args)
  # Discard patients that were filtered in processing
  study_data_complete <-  study_data_complete %>% rownames_to_column("patient_id") %>% 
    filter(patient_id %in% colnames(betas_complete$cpg_beta$cpg_beta_df))
  rownames(study_data_complete) <- study_data_complete$patient_id
 
  
  df_cpg <- data.frame(cpg_id = tdhia::manifest_v1A2_design_scores$cpg_id) %>% filter(!is.na(cpg_id))
  df_icr <- data.frame(icr_id = tdhia::manifest_v1A2_design_scores$icr_id) %>% filter(!is.na(icr_id))
  
  # Parse model string, test if response is binomial or not
  if (is.null(family)) {
    out <- tdhia::parse_lm_formula(model_str)
    if (length(unique(study_data[[out$response]])) ==2 ) {
      family <- "binomial"
    } else {family <- "gaussian"}
    verbosecat(sprintf("Inferred model from response var: %s\n", family))
  }
  # study_data %>% select(out$covariates, out$response)
  # df_cpg_glm %>% select()
  
  
  df_cpg_glm_path = file.path(data_cache_path, "df_cpg_glm.rds")
  if (!file.exists(df_cpg_glm_path) || overwrite_cache) {
  df_cpg_glm <- tdhia::imprintome_glm(
    model_str = model_str, study_data = study_data,
    betas = betas_na$cpg_beta$cpg_beta_df, family = family, m_value_transform = m_value_transform,
    n_p_adj = nrow(betas_na$icr_beta$icr_beta_df), db_flag = T, rm.na.all = !impute_na, 
    verbose = verbose, impute_na = impute_na, max_p_val = 0.05, n.cores = n.cores)
  saveRDS(object = df_cpg_glm,  file = df_cpg_glm_path)
  } else { 
    df_cpg_glm <- readRDS(file = df_cpg_glm_path)
  }
  # Add results to master cpg dataframe
  df_cpg <- df_cpg %>% left_join( 
    df_cpg_glm$imp_site %>% select(Variable, Estimate, Statistic, P_VAL, ADJ_P_VAL, Family) %>%
      rename(cpg_id = Variable, cpg_glm_estimate = Estimate, cpg_glm_statistic = Statistic,
             cpg_glm_raw_pval = P_VAL, cpg_glm_adj_pval = ADJ_P_VAL, cpg_glm_family = Family), by = join_by(cpg_id)) 
 
  
  
  # df_icr_glm_path = file.path(data_cache_path, "df_icr_glm.rds")
  # if (!file.exists(df_icr_glm_path)) {
  # df_icr_glm <- tdhia::imprintome_glm(
  #   model_str = model_str, study_data = study_data,
  #   betas = betas_na$icr_beta$icr_beta_df, family = family, m_value_transform = m_value_transform,
  #   n_p_adj = nrow(betas_na$icr_beta$icr_beta_df), db_flag = T, rm.na.all = !impute_na, 
  #   verbose = verbose, impute_na = impute_na, max_p_val = 0.05, n.cores = n.cores)
  # saveRDS(object = df_icr_glm,  file = df_icr_glm_path)
  # } else { 
  #   df_icr_glm <- readRDS(file = df_icr_glm_path)
  # }
  
  df_icr_skat_path = file.path(data_cache_path, "df_icr_skat.rds")
  if (!file.exists(df_icr_skat_path) || overwrite_cache) {
   df_icr_skat <- tdhia::skat_icr_test(
     cpg_betas = as.data.frame(betas_complete$cpg_beta$cpg_beta_df),
     df_study = study_data_complete,  response = parsed_model$response,
     predictors = parsed_model$covariates,  method = "optimal.adj",
     out_type = "C",  icr_ids = NULL,  min_cpg = 3,  db_flag = FALSE,
     m_value_transform = m_value_transform,  scaling = TRUE,  verbose = verbose,
     n.cores = n.cores)
   saveRDS(object = df_icr_skat,  file = df_icr_skat_path)
  } else { 
    df_icr_skat <- readRDS(file = df_icr_skat_path)
  }
  df_icr <- df_icr %>% left_join( 
    df_icr_skat %>% rename(skat_n_cpg = n_cpg), by = join_by(icr_id)) 
  
  
 
     
  
  df_icr_pcr_path = file.path(data_cache_path, "df_icr_pcr.rds")
  if (!file.exists(df_icr_pcr_path) || overwrite_cache) {
   df_icr_pcr <- tdhia::pc_regression_test(
    cpg_beta = Matrix::t(betas_complete$cpg_beta$cpg_beta_df) %>% as.data.frame(),
    m_value_transform = m_value_transform,  data_norm_type = "n1", pct_variance = 0.80,
    df_study = study_data_complete %>% rename(Patient_ID = patient_id),
    outcome = parsed_model$response,  covariates = parsed_model$covariates,
    Patient_ID = "Patient_ID",  family = "gaussian",  icr_ids = NULL,
    min_cpg = 3,  verbose = verbose,  n.cores = 1)
   saveRDS(object = df_icr_pcr,  file = df_icr_pcr_path)
  } else { 
    df_icr_pcr <- readRDS(file = df_icr_pcr_path)
  }
  df_icr <- df_icr %>% left_join( 
    df_icr_pcr %>% rename(icr_id = ICR_id) %>% 
      rename_with(~ paste0("pcr_", .), -icr_id), by = join_by(icr_id)) 
  
 
  
  df_cpg_limma_path = file.path(data_cache_path, "df_cpg_limma.rds")
  if (!file.exists(df_cpg_limma_path) || overwrite_cache) {
   df_cpg_limma <- tdhia::cpg_dml_test(
     df_study = study_data_complete, predictors = c(out$response,out$covariatesr), 
     cpg_beta = betas_complete$cpg_beta$cpg_beta_df,
     pvalue_threshold = 0.001, db_flag = T, sample_name = "patient_id", correlation_check = F,
     m_value_transform = m_value_transform,
     beadchip_correction = F, verbose = T)
   saveRDS(object = df_cpg_limma,  file = df_cpg_limma_path)
  } else { 
    df_cpg_limma <- readRDS(file = df_cpg_limma_path)
  }
  df_cpg <- df_cpg %>% left_join( 
    df_cpg_limma$df_dml %>% select(CpG_Probe, logFC, AveExpr, P.Value, t, P.Value, adj.P.Val  ) %>%
      rename(cpg_id = CpG_Probe, cpg_limma_logfc = logFC, cpg_limma_avg_expr = AveExpr,
             cpg_limma_raw_pvalue = P.Value, cpg_limma_t = t, cpg_limma_adj_pvalue = adj.P.Val), by = join_by(cpg_id)) 
  
  
  
  df_icr_lancaster_path = file.path(data_cache_path, "df_icr_lancaster.rds")
  if (!file.exists(df_icr_lancaster_path) || overwrite_cache) {
   df_icr_lancaster <-  tdhia::icr_dmr_test(
     df_dml = df_cpg_limma$df_dml, chr_lens = df_cpg_limma$chr_lens,
     pval_threshold = 0.05, fdr_sig_threshold = 0.0001, verbose = T, db_flag = F)
   saveRDS(object = df_icr_lancaster,  file = df_icr_lancaster_path)
  } else { 
    df_icr_lancaster <- readRDS(file = df_icr_lancaster_path)
  }
  df_icr <- df_icr %>% left_join( 
    df_icr_lancaster$ICR_summary %>% rename(icr_id = ICR_id) %>% 
      rename_with(~ paste0("pcr_", .), -icr_id), by = join_by(icr_id)) 
  
   
   
}
