



tdha_stat_tests <- function(
    model_str, beta_data, study_data, transform_to_m = T, cache_path, 
    model_prefix = "", impute_na = TRUE, n.cores = max(c(parallel::detectCores()-1, 1))) {
  
  
  # Loop through each model
  # Process all stats tests
  # 1) cpg_glm
  # 2) icr_glm
  # 2) icr_skat
  # 3) icr_pcr
  # 4) cpg_limma
  # 5) icr_lancaster
  
  if (impute_na) {
    beta_data_na <- 
  } else {beta_data_na = beta_data}
 
      
      
  df_cpg_glm <- tdhia::imprintome_glm(
    model_str = model_str, study_data = study_data,
    betas = beta_data$cpg_beta$cpg_beta_df, family = "binomial", m_value_transform = m_value_transform,
    n_p_adj = nrow(data$icr_beta$icr_beta_df), db_flag = FALSE, rm.na.all = !impute_na, 
    verbose = TRUE, impute_na = impute_na, max_p_val = 0.05, n.cores = n.cores)
  
  
  df_icr_glm <- tdhia::imprintome_glm(
    model_str = model_str, study_data = study_data,
    betas = beta_data$icr_beta$icr_beta_df, family = "binomial", m_value_transform = m_value_transform,
    n_p_adj = nrow(data$icr_beta$icr_beta_df), db_flag = FALSE, rm.na.all = !impute_na, 
    verbose = TRUE, impute_na = impute_na, max_p_val = 0.05, n.cores = n.cores)
  
  
  
  
  
  
   df_icr_skat <- tdhia::skat_icr_test(
     cpg_betas = as.data.frame(data_no_na$cpg_beta$cpg_beta_df),
     df_study = study_metadata,  response = "disease_state",
     predictors = c("Sex", "Age", "Ethnicity"),  method = "optimal.adj",
     out_type = "D",  icr_ids = NULL,  min_cpg = 3,  db_flag = FALSE,
     m_value_transform = m_value_transform,  scaling = TRUE,  verbose = TRUE,
     n.cores = n.cores)
     
     
   df_icr_pcr <- tdhia::pc_regression_test(
    cpg_beta = Matrix::t(data_no_na$cpg_beta$cpg_beta_df) %>% as.data.frame(),
    m_value_transform = m_value_transform,  data_norm_type = "n1", pct_variance = 0.80, n_pcs = 1,
    df_study = study_metadata %>% rename(Patient_ID = Patient.ID),
    outcome = "disease_state",  covariates = c("Sex"),
    Patient_ID = "Patient_ID",  family = "binomial",  icr_ids = NULL,
    min_cpg = 3,  verbose = TRUE,  n.cores = 1)
  
   
   df_cpg_limma <- tdhia::cpg_dml_test(
     df_study = study_metadata, predictors = c("disease_state", "Age", "Sex", "Ethnicity"), cpg_beta = data$cpg_beta$cpg_beta_df,
     pvalue_threshold = 0.001, db_flag = F, sample_name = "Patient.ID", correlation_check = F,
     m_value_transform = m_value_transform,
     beadchip_correction = F, verbose = T, write_plots = F,
     output_dir_path = getwd())
   
   df_icr_lancaster <-  tdhia::icr_dmr_test(
     df_dml = cpg_dml_out$df_dml, chr_lens = cpg_dml_out$chr_lens,
     pval_threshold = 0.05, fdr_sig_threshold = 0.0001, verbose = T, db_flag = F)
   
   
   
   
}
