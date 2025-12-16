


#' skat_icr_test
#'
#' performs skat test for each icr test.
#'
#' @param cpg_betas cpg beta matrix (cpg_id as rows x sample_id as columns)
#' @param df_study dataframe of sample variables to be used in linear models, nrows = sample size.
#' @param response string, column name of response variable located in df_study.
#' @param predictors vector list of column names of predictors for linear model.
#' @param icr_ids vector list of icr_ids to be tested. Default = NULL, tests all 
#' icrs that are covered by the input cpg beta matrix.
#' @param apply_pca boolean, when true, applies skat to the first ncomps of a
#' PCA of the cpg_betas matrix. Default=FALSE.
#' @param ncomp numebr of PCA components to use when apply_pca=TRUE.
#' @returns test
#' @export
skat_icr_test <- function(cpg_betas, df_study, response, predictors,  icr_ids = NULL, 
                          apply_pca = F, ncomp = 10, db_flag = T, transform_m_values = T){
  if(db_flag) save(list = ls(all.names = TRUE), file = "skat_icr_test.RData")
  # load(file = "skat_icr_test.RData")
  
  cpg_mapping <- tdhia::manifest_v1A2_design_scores %>% 
    dplyr::select(cpg_id, icr_id)
  
  # if icr_ids no supplied, scan for all icr_ids covered with cpg_ids
  if (is.null(icr_ids)) {
    icr_ids <- cpg_mapping %>% 
      dplyr::filter(cpg_id %in% rownames(cpg_betas)) %>% dplyr::pull(icr_id) %>%
      unname() %>% unique()
  }
  
  # Remove all samples with NA values from study data (not supported with SKAT)
  df_study <- df_study %>% select(all_of(c(response, predictors))) %>% tidyr::drop_na()
  
  # Make the samples match and order them the same
  shared_sample_ids <- intersect(colnames(cpg_betas),  rownames(df_study))
  cpg_betas <- cpg_betas[,shared_sample_ids]
  df_study <- df_study[rownames(df_study) %in% shared_sample_ids, ]
  
  # Order df_study the same sample order as cpg_beta
  df_study <- df_study[match(x = colnames(cpg_betas), table = rownames(df_study)),]
  
  # Dataframe to store results of test
  df_results <- data.frame(icr_id = icr_ids, skat_p_value = NA, n_cpg = NA)
  
  model_str = paste0(response , " ~ ", paste(predictors, collapse = " + "))
  
  # parfun = function(x) {
  #   skat_parfun(icr_id_ = df_results$icr_id[x], cpg_betas, df_study, model_str, cpg_mapping)
  # }

  
  for (n in 1:nrow(df_results)) {
    # print(n)
    # Get list of cpgs for a given ICR
    subset_cpg_ids <- cpg_mapping %>% filter(icr_id == df_results$icr_id[n]) %>% pull(cpg_id)
    
    # Subset cpg beta matrix to only those contained within ICR
    tZ <- cpg_betas %>% rownames_to_column("cpg_id") %>% filter(cpg_id %in% subset_cpg_ids) %>% 
      column_to_rownames("cpg_id")
    # 
    if (transform_m_values) tZ = sesame::BetaValueToMValue(tZ)
    
    # Construct model string
    model_str = paste0(response , " ~ ", paste(predictors, collapse = " + "))
    # Calculate SKAT null model
    skat_null <- SKAT::SKAT_Null_Model(formula =  as.formula(model_str), data = df_study, out_type="C",
                           n.Resampling = 0, Adjustment = TRUE)
    # SKAT observed model
    skat_out<-SKAT::SKAT(Z = Matrix::t(tZ), obj = skat_null, kernel = "linear")
    
    df_results$skat_p_value[n] <- skat_out$p.value
    df_results$n_cpg[n] = nrow(tZ)
  }
  
  # Calculate adjusted p-value and q-value
  df_results$adj_p_value <- p.adjust(p = df_results$skat_p_value, method = "fdr")
  df_results$q_value <- qvalue::qvalue(p = df_results$skat_p_value, fdr.level = 0.05)$qvalues
  
  
  return(df_results)
  
}


skat_parfun <- function(icr_id_, cpg_betas, df_study, model_str, cpg_mapping) {
  # Get list of cpgs for a given ICR
  subset_cpg_ids <- cpg_mapping %>% filter(icr_id == icr_id_) %>% pull(cpg_id)
  
  # Subset cpg beta matrix to only those contained within ICR
  tZ <- cpg_betas %>% rownames_to_column("cpg_id") %>% filter(cpg_id %in% subset_cpg_ids) %>% 
    column_to_rownames("cpg_id")
  
  # Construct model string
  
  # Calculate SKAT null model
  skat_null <- SKAT::SKAT_Null_Model(formula =  as.formula(model_str), data = df_study, out_type="D",
                                     n.Resampling = 0, Adjustment = TRUE)
  # SKAT observed model
  skat_out<-SKAT::SKAT(Z = Matrix::t(tZ), obj = skat_null, kernel = "linear")
  
  out <- data.frame(icr_id = icr_id_, skat_p_value = skat_out$p.value, n_cpg = nrow(tZ))
  
  return(out)
}
