#' Run CpG and ICR Methylation Association Tests
#'
#' Runs generalized linear models (GLMs) and limma tests at CpG level, and
#' GLMs, SKAT, principal component regression, and Lancaster combination tests
#' at imprinting control region (ICR) level. Combines results into annotated
#' tables and counts significant results for each adjusted p-value column.
#'
#' @param primary_var Character string naming a column in `study_data`. Used
#'   as the response in GLM, SKAT, and principal component regression, and as
#'   the first predictor in the limma analysis.
#' @param predictor_vars Nonempty character vector of additional column names
#'   in `study_data`. All are included in GLM and limma models. The current
#'   implementation uses only `predictor_vars[-1]` as covariates for SKAT and
#'   principal component regression.
#' @param betas List returned by [tdhia_pipeline()], including `input_args`
#'   and `probe_beta`. The pipeline is rerun from these inputs to prepare
#'   methylation data for the tests; this argument is not a beta matrix.
#' @param study_data Data frame with one row per sample and columns named by
#'   `primary_var` and `predictor_vars`. Row names must be sample identifiers
#'   matching the methylation matrix column names and pipeline sample IDs.
#' @param m_value_transform Logical; whether to transform beta values to
#'   M-values in the underlying tests. Defaults to `TRUE`.
#' @param data_cache_path Character string giving an existing, writable
#'   directory for intermediate RDS results. The directory is not created.
#' @param model_prefix Character string. Currently unused, including in cache
#'   file names.
#' @param impute_na Logical; whether the CpG and ICR GLMs impute missing values
#'   using MICE. If `FALSE`, missing values are handled by removal in the GLM
#'   helper. This option does not control the other tests.
#' @param n.cores Number of cores passed to GLM and SKAT helpers. Defaults to
#'   the detected core count minus four, with a minimum of one. Principal
#'   component regression is always run with one core.
#' @param verbose Logical; whether to print progress from this function and
#'   helpers that receive this argument. Limma and Lancaster helpers are
#'   currently called with verbose output enabled regardless of this value.
#' @param family Character string specifying the GLM family, such as
#'   `"gaussian"` or `"binomial"`, also passed to principal component regression.
#'   If `NULL`, uses `"binomial"` when `study_data[[primary_var]]` has exactly
#'   two unique values, and `"gaussian"` otherwise. Missing values count as a
#'   unique value in this inference. SKAT is always called with `out_type = "C"`
#'   (continuous response), independently of this argument.
#' @param db_flag Logical; whether to save the function's initial environment
#'   to `tdhia_stat_tests.RData` in the working directory. Defaults to `FALSE`.
#'   Pipeline and GLM helpers are called with their own debug flags enabled
#'   even when this argument is `FALSE`.
#' @param overwrite_stats_cache Logical; whether to recompute cached CpG GLM, SKAT,
#'   principal component regression, limma, and Lancaster results. Defaults
#'   to `FALSE`. Existing ICR GLM results are always reused, even when `TRUE`.
#' @param add_genomic_metadata Logical; whether to join genomic annotations
#'   from [add_metadata_to_imp_sites()] to the CpG and ICR result tables.
#'   Defaults to `TRUE`.
#'
#' @details
#' The GLM formula is `primary_var ~ beta + predictor_vars`, with the named
#' predictors joined by `+` and `beta` representing each methylation site.
#' GLMs use pipeline output with failed beta measurements set to missing.
#' The other tests use samples with complete study variables and pipeline
#' output that retains failed beta measurements. Samples removed by pipeline
#' filtering are excluded from the corresponding prepared study data.
#'
#' SKAT uses `method = "optimal.adj"`, scaling, and at least three CpGs per
#' ICR. Principal component regression also requires at least three CpGs and
#' retains components explaining 80 percent of variance. Both GLM analyses
#' use the number of ICRs for the number of p-value adjustment comparisons.
#'
#' Intermediate test results are cached using file names derived only from
#' `primary_var` and the test name. Changes to data, predictors, or other
#' settings do not invalidate caches. Use a separate cache directory for
#' each analysis configuration. Pipeline preparation runs even when test
#' results are read from cache.
#'
#' @return A named list of four data frames:
#' \describe{
#'   \item{df_cpg}{CpGs assigned to ICRs in `manifest_v1A2_design_scores`,
#'     identified by `primary_var` and `cpg_id`, with `cpg_glm_*` and
#'     `cpg_limma_*` statistics and optional genomic metadata.}
#'   \item{df_cpg_summary}{One row containing `primary_var` and counts of
#'     values below 0.05 in each numeric adjusted p-value column of `df_cpg`,
#'     ignoring missing values.}
#'   \item{df_icr}{ICRs in the manifest, identified by `primary_var` and
#'     `icr_id`, with `icr_glm_*`, `icr_skat_*`, `icr_pcr_*`, and `icr_lanc_*`
#'     results, CpG-level summaries, and optional genomic metadata. CpG
#'     summaries include minimum p-values and signed values of greatest
#'     absolute magnitude for effect estimates and test statistics.}
#'   \item{df_icr_summary}{One row containing `primary_var` and counts of
#'     values below 0.05 in each numeric adjusted p-value column of `df_icr`,
#'     ignoring missing values.}
#' }
#' Results are left-joined to manifest identifiers, so untested sites can
#' have missing statistics.
#'
#' @seealso [tdhia_pipeline()], [imprintome_glm()], [skat_icr_test()],
#'   [pc_regression_test()], [cpg_dml_test()], [icr_dmr_test()]
#'
#' @export
tdhia_stat_tests <- function( 
    primary_var, predictor_vars, betas, study_data, m_value_transform = T, 
    data_cache_path, model_prefix = "", impute_na = TRUE, 
    n.cores = max(c(parallel::detectCores()-4, 1)), verbose = T, 
    family = NULL,db_flag = FALSE, overwrite_stats_cache = F,
    overwrite_betas_cache = F, add_genomic_metadata = T) {
  
  if (db_flag) {save(list = ls(all.names = TRUE), file = "tdhia_stat_tests.RData")}
  # load(file = "tdhia_stat_tests.RData")
  verbosecat <- \(x) if (verbose) cat(x)
  dir.create(data_cache_path, showWarnings = F,recursive = T)
  
  # Evaluate the analysis only when its cached result cannot be reused.
  cache_result <- function(path, expr, overwrite = FALSE) {
    if (file.exists(path) && !overwrite) return(readRDS(path))
    result <- force(expr)
    saveRDS(result, path)
    result
  }
  
  # Set model variables
  model_str <- paste0(primary_var, " ~ beta + ", paste0(predictor_vars,collapse = " + "))
  
  # Extract columns used in study data
  study_columns <- c(unname(primary_var), predictor_vars)
  study_columns <- study_columns[study_columns != "beta"]
  
  beta_path_na <-  file.path(data_cache_path, paste0(model_prefix, "_", primary_var, "_na_betas_and_study_data.rds"))
  beta_path_complete <-  file.path(data_cache_path, paste0(model_prefix, "_", primary_var, "_complete_betas_and_study_data.rds"))
  
  
  # Betas and study data with NAs                                          ####
  #_____________________________________________________________________________
  if (!file.exists(beta_path_na) || overwrite_betas_cache) {
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
    # Package for export
    data_na <- list(betas = betas_na, study_data = study_data_na)
    saveRDS(data_na, file = beta_path_na)
    remove(betas_na, study_data_na)
    
  } else { data_na <- readRDS(beta_path_na) }
  
  
  
  # Betas and study data without NAs                                        ####
  #_____________________________________________________________________________
  if (!file.exists(beta_path_complete) || overwrite_betas_cache) {
    study_data_complete <- study_data %>% select(all_of(study_columns)) %>% na.omit()
    betas_complete <- betas
    # Only include patient columns that are found in complete records in study data
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
    # Package for export
    data_complete = list(betas = betas_complete, study_data = study_data_complete)
    saveRDS(data_complete, file = beta_path_complete)
    remove(study_data_complete, betas_complete)
    
  } else { data_complete <- readRDS(beta_path_complete) }
  
  
  # Initialize result grame for cpg and icr level analysis
  cpg_ids = tdhia::manifest_v1A2_design_scores %>% filter(!is.na(icr_id)) %>% 
    pull(cpg_id) %>% unique()
  df_cpg <- data.frame(primary_var = rep(primary_var, length(cpg_ids)), cpg_id = cpg_ids)
  icr_ids = tdhia::manifest_v1A2_design_scores  %>% filter(!is.na(icr_id)) %>%
    pull(icr_id) %>% unique()
  df_icr <- data.frame(primary_var = rep(primary_var, length(icr_ids)), icr_id = icr_ids)
  rm(cpg_ids, icr_ids)
  
  # Parse model string, test if response is binomial or not
  if (is.null(family)) {
    if (length(unique(study_data[[primary_var]])) ==2 ) {
      family <- "binomial"
    } else {family <- "gaussian"}
    verbosecat(sprintf("Inferred model from response var: %s\n", family))
  }
  
  
  # Helper function for summarizing cpg level results in icr table
  # Return max abs value with sign preserved
  max_mag_sign = \(x) sign(x[which.max(abs(x))[1]]) * max(abs(x))
  
  
  # CpG GLM                                                            ##########  
  #_____________________________________________________________________________
  df_cpg_glm_path = file.path(data_cache_path, paste0(model_prefix, "_", primary_var, "_df_cpg_glm.rds"))
  df_cpg_glm <- cache_result(df_cpg_glm_path, {
    tdhia::imprintome_glm(
      model_str = model_str, study_data = data_na$study_data,
      betas = data_na$betas$cpg_beta$cpg_beta_df, family = family, m_value_transform = m_value_transform,
      n_p_adj = nrow(data_na$betas_na$icr_beta$icr_beta_df), db_flag = T, rm.na.all = !impute_na,
      verbose = verbose, impute_na = impute_na, max_p_val = 0.05, n.cores = n.cores)
  }, overwrite = overwrite_stats_cache)
  # Add results to master cpg dataframe
  df_cpg_glm_formatted <- df_cpg_glm$imp_site %>% 
    select(Variable, Estimate, Statistic, P_VAL, ADJ_P_VAL, Family) %>%
    rename(cpg_id = Variable, cpg_glm_estimate = Estimate, cpg_glm_statistic = Statistic,
           cpg_glm_raw_pval = P_VAL, cpg_glm_adj_pval = ADJ_P_VAL, cpg_glm_family = Family) %>%
    left_join(  tdhia::mapping_cpg_icr_ids %>% select(ICR_id, CpG_id) %>% distinct() %>% 
                  rename(icr_id = ICR_id, cpg_id = CpG_id), by = join_by(cpg_id))
  df_cpg <- df_cpg %>% left_join( df_cpg_glm_formatted %>% select(-icr_id), by = join_by(cpg_id)) 
  # Add cpg_glm results to icr table
  df_icr <- df_icr %>% left_join(
    df_cpg_glm_formatted %>% group_by(icr_id) %>% summarize(
      cpg_glm_estimate_max = max_mag_sign(cpg_glm_estimate), 
      cpg_glm_statistic_max = max_mag_sign(cpg_glm_statistic),
      cpg_glm_raw_pval_min = min(cpg_glm_raw_pval), cpg_glm_adj_pval_min = min(cpg_glm_adj_pval),
      cpg_glm_family = cpg_glm_family[1]
    ), by = join_by(icr_id) )
  
  
  # ICR GLM                                                    #################
  #_____________________________________________________________________________
  df_icr_glm_path = file.path(data_cache_path, paste0(model_prefix, "_", primary_var, "df_icr_glm.rds"))
  df_icr_glm <- cache_result(df_icr_glm_path, {
    tdhia::imprintome_glm(
      model_str = model_str, study_data = data_na$study_data,
      betas = data_na$betas$icr_beta$icr_beta_df, family = family, m_value_transform = m_value_transform,
      n_p_adj = nrow(data_na$betas$icr_beta$icr_beta_df), db_flag = T, rm.na.all = !impute_na,
      verbose = verbose, impute_na = impute_na, max_p_val = 0.05, n.cores = n.cores)
  }, overwrite = FALSE)
  # Add cpg level results
  df_icr_glm_formatted <- df_icr_glm$imp_site %>% 
    select(Variable, Estimate, Statistic, P_VAL, ADJ_P_VAL, Family) %>%
    rename(icr_id = Variable, icr_glm_estimate = Estimate, icr_glm_statistic = Statistic,
           icr_glm_raw_pval = P_VAL, icr_glm_adj_pval = ADJ_P_VAL, icr_glm_family = Family)
  # Add icr level results
  df_icr <- df_icr %>% left_join(df_icr_glm_formatted, by = join_by(icr_id) )
  
  
  
  
  # ICR SKAT                                                         ##########  
  #_____________________________________________________________________________
  df_icr_skat_path = file.path(data_cache_path, paste0(model_prefix, "_", primary_var, "df_icr_skat.rds"))
  df_icr_skat <- cache_result(df_icr_skat_path, {
    tdhia::skat_icr_test(
      cpg_betas = as.data.frame(data_complete$betas$cpg_beta$cpg_beta_df),
      df_study = data_complete$study_data,  response = unname(primary_var),
      predictors = predictor_vars[-1],  method = "optimal.adj",
      out_type = "C",  icr_ids = NULL,  min_cpg = 3,  db_flag = FALSE,
      m_value_transform = m_value_transform,  scaling = TRUE,  verbose = verbose,
      n.cores = n.cores)
  }, overwrite = overwrite_stats_cache)
  # Add ICR level results
  df_icr <- df_icr %>% left_join( 
    df_icr_skat %>% rename(skat_n_cpg = n_cpg) %>% rename_with(~ paste0("icr_", .), -icr_id), 
    by = join_by(icr_id)) 
  
  
  
  
  # ICR PCA Regression                                               ##########  
  #_____________________________________________________________________________
  df_icr_pcr_path = file.path(data_cache_path, paste0(model_prefix, "_", primary_var, "df_icr_pcr.rds"))
  df_icr_pcr <- cache_result(df_icr_pcr_path, {
    tdhia::pc_regression_test(
      cpg_beta = Matrix::t(data_complete$betas$cpg_beta$cpg_beta_df) %>% as.data.frame(),
      m_value_transform = m_value_transform,  data_norm_type = "n1", pct_variance = 0.80,
      df_study = data_complete$study_data %>% rename(Patient_ID = patient_id),
      outcome = unname(primary_var),  covariates = predictor_vars[-1],
      Patient_ID = "Patient_ID",  family = family,  icr_ids = NULL,
      min_cpg = 3,  verbose = verbose,  n.cores = 1)
  }, overwrite = overwrite_stats_cache)
  # Add ICR level results
  df_icr <- df_icr %>% left_join( 
    df_icr_pcr %>% rename(icr_id = ICR_id, adj_pval = adj_p_value ) %>% 
      rename_with(~ paste0("icr_pcr_", .), -icr_id), by = join_by(icr_id)) 
  
  
  # CpG Limma                                                         ##########  
  #_____________________________________________________________________________
  df_cpg_limma_path = file.path(data_cache_path, paste0(model_prefix, "_", primary_var, "df_cpg_limma.rds"))
  df_cpg_limma <- cache_result(df_cpg_limma_path, {
    tdhia::cpg_dml_test(
      df_study = data_complete$study_data, predictors = c(primary_var, predictor_vars),
      cpg_beta = data_complete$betas$cpg_beta$cpg_beta_df,
      pvalue_threshold = 0.001, db_flag = F, sample_name = "patient_id", correlation_check = F,
      m_value_transform = m_value_transform,
      beadchip_correction = F, verbose = T)
  }, overwrite = overwrite_stats_cache)
  # Add cpg level results
  df_cpg_limma_formatted <- df_cpg_limma$df_dml %>% 
    select(CpG_Probe, logFC, AveExpr, P.Value, t, P.Value, adj.P.Val  ) %>%
    rename(cpg_id = CpG_Probe, cpg_limma_logfc = logFC, cpg_limma_avg_expr = AveExpr,
           cpg_limma_raw_pval = P.Value, cpg_limma_t = t, cpg_limma_adj_pval = adj.P.Val) %>% 
    left_join(  tdhia::mapping_cpg_icr_ids %>% select(ICR_id, CpG_id) %>% distinct() %>% 
                  rename(icr_id = ICR_id, cpg_id = CpG_id), by = join_by(cpg_id))
  df_cpg <- df_cpg %>% left_join( select(df_cpg_limma_formatted,-icr_id), by = join_by(cpg_id)) 
  # Add cpg limma results to icr results
  df_icr <- df_icr %>% left_join(
    df_cpg_limma_formatted %>% group_by(icr_id) %>% summarize(
      cpg_limma_logfc_max  = max_mag_sign(cpg_limma_logfc ), 
      cpg_limma_avg_expr_max  = max_mag_sign(cpg_limma_avg_expr ),
      cpg_limma_t_max  = max_mag_sign(cpg_limma_t),
      cpg_limma_raw_pval_min = min(cpg_limma_raw_pval ), 
      cpg_limma_adj_pval_min = min(cpg_limma_adj_pval)), by = join_by(icr_id))
  
  
  # ICR lancaster                                                     ##########  
  #_____________________________________________________________________________
  df_icr_lancaster_path = file.path(data_cache_path, paste0(model_prefix, "_", primary_var, "df_icr_lancaster.rds"))
  df_icr_lancaster <- cache_result(df_icr_lancaster_path, {
    tdhia::icr_dmr_test(
      df_dml = df_cpg_limma$df_dml, chr_lens = df_cpg_limma$chr_lens,
      pval_threshold = 0.05, fdr_sig_threshold = 0.0001, verbose = T, db_flag = F)
  }, overwrite = overwrite_stats_cache)
  df_icr <- df_icr %>% left_join( 
    df_icr_lancaster$ICR_summary %>% rename(icr_id = ICR_id, combined_adj_pval = FDR) %>% 
      rename_with(~ paste0("icr_lanc_", .), -icr_id), by = join_by(icr_id)) 
  #  cpg level results
  df_cpg_summary <- df_cpg %>% select(contains("adj_pval")) %>%
    summarise(across(where(is.numeric), ~ sum(.x < 0.05, na.rm = TRUE))) %>% 
    mutate(primary_var = primary_var, .before = 1)
  #  icr level results
  df_icr_summary <- df_icr %>% select(contains("adj_pval")) %>%
    summarise(across(where(is.numeric), ~ sum(.x < 0.05, na.rm = TRUE))) %>%
    mutate(primary_var = primary_var, .before = 1)
  
  
  # Add metadata to results tables
  if (add_genomic_metadata) {
    df_cpg <-  df_cpg %>% left_join(tdhia::add_metadata_to_imp_sites(
      df_cpg$cpg_id, imp_type = "cpg"), by = join_by("cpg_id"), keep = F,multiple = "first")
    df_icr <-  df_icr %>% left_join(tdhia::add_metadata_to_imp_sites(
      df_icr$icr_id, imp_type = "icr"), by = join_by("icr_id"), keep = F,multiple = "first")
  }
  
  
  # Record column prefixes for each analysis
  prefix = list()
  prefix$cpg_glm = "cpg_glm"
  prefix$df_icr_glm = "icr_glm"
  prefix$df_icr_skat = "icr_skat"
  prefix$df_icr_pcr = "icr_pcr"
  prefix$df_cpg_limma = "cpg_limma"
  prefix$df_icr_lancaster = "icr_lancaster"
  
  
  
  
  
  
  # Export
  out <- list(df_cpg = df_cpg, df_cpg_summary = df_cpg_summary,
              df_icr = df_icr, df_icr_summary = df_icr_summary)
  
  # For each list in out, apply attributes
  for (n in seq_along(out)) {
    attr(out[[n]], "prefix") <- prefix
    attr(out[[n]], "labels") <- tdhia_stat_pretty_labels(out[[n]])
    
  }
  
  return(out)
}


# returns colnames and pretty labels
tdhia_stat_pretty_labels <- function(df) {
  
  col_names = colnames(df)
  pretty_labels = col_names
  
  pretty_labels <- str_replace_all(pretty_labels, "_", " ")
  
  pretty_labels <- str_replace_all(pretty_labels, "pval", "p-Value")
  pretty_labels <- stringr::str_to_title(pretty_labels)
  
  names(col_names) <- pretty_labels
  
  return(col_names)
}