

#' skat_icr_test
#'
#' Performs kernel-regression-based association tests for methylation state for
#' a collection of cpg sites. Both continuous and dichotomous response variables
#' are supported.
#'
#' @param cpg_betas cpg beta matrix (cpg_id as rows x sample_id as columns).
#'  Values are assumed to be on beta scale (0-1). Use m_value_transform to
#'   convert to m-values (recommended).
#' @param df_study dataframe of sample associated data to be used in linear
#' models, (nrows = sample size). Columns should include those specified with
#' 'response' and 'predictors' input arguments.
#' @param response string, column name of response variable located in df_study.
#' @param predictors vector of strings of column names of predictors that are
#' located in df_study.
#' @param method
#'   davies: exact p-value for skat method (heterogeneous)
#'   optimal.adj: skat-O unified approach, combination of SKAT and burden test
#'    (default).
#' @param out_type character, specifies type of variable for response column.
#' "C": continuous (default).
#' "D":  dichotomous.
#' @param icr_ids vector of strings of icr_ids to be tested. Default = NULL,
#' tests all icrs that are covered by the input cpg beta matrix.
#' @param min_cpg minimum number of cpg sites for an icr to be included in results.
#' @param m_value_transform boolean, when true, transform beta values to
#' m-values to control heteroskedasticity. Default = TRUE.
#' @param verbose boolean, when TRUE, prints progress to command line (default = TRUE).
#' @param n.cores integer, number of cores for processing (default = 1).
#' @returns test
#' @export
#' @author Kate Everly, Bruce Corliss
skat_icr_test <- function(cpg_betas, df_study, response, predictors,
                          method = "optimal.adj",
                          out_type="C", icr_ids = NULL,
                          min_cpg = 3, db_flag = T, m_value_transform = T,
                          scaling = T,
                          verbose = T, n.cores = 1){
  if(db_flag) save(list = ls(all.names = TRUE), file = "skat_icr_test.RData")
  # load(file = "skat_icr_test.RData")
  verbosecat <- function(x) if(verbose) cat(x)

  cpg_mapping <- tdhia::manifest_v1A2_design_scores %>%
    dplyr::select(cpg_id, icr_id)

  # if icr_ids no supplied, scan for all icr_ids covered with cpg_ids
  if (is.null(icr_ids)) {
    icr_ids <- cpg_mapping %>%
      dplyr::filter(cpg_id %in% rownames(cpg_betas)) %>% dplyr::pull(icr_id) %>%
      unname() %>% unique()
  }

  verbosecat("> Remove rows in df_study that contain NA...\n")
  # Remove all samples with NA values from study data (not supported with SKAT)
  df_study <- df_study %>% dplyr::select(all_of(c(response, predictors))) %>% tidyr::drop_na()

  # Make the samples match and order them the same
  verbosecat("> Forcing sample id order to match between data and df_study.\n")
  shared_sample_ids <- intersect(colnames(cpg_betas),  rownames(df_study))
  cpg_betas <- cpg_betas[,shared_sample_ids]
  df_study <- df_study[rownames(df_study) %in% shared_sample_ids, ]

  # Order df_study the same sample order as cpg_beta
  df_study <- df_study[match(x = colnames(cpg_betas), table = rownames(df_study)),]

  # Dataframe to store results of test
  df_results <- data.frame(icr_id = icr_ids, skat_p_value = NA, n_cpg = NA)

  model_str = paste0(response , " ~ ", paste(predictors, collapse = " + "))
  verbosecat(sprintf("> Model: %s\n", model_str))

  # packaged skat into function for single ICR to make parallel conversion easier in future
  if (n.cores==1){ verbosecat("> SKAT processing on single core.\n")

    out = list()
    for (n in 1:length(icr_ids)) {

      out[[n]] <- skat_single_icr(icr_id = icr_ids[n], cpg_betas, df_study,
                                  model_str, cpg_mapping, m_value_transform, scaling = scaling, method = method, out_type = out_type)
    }; df_results = do.call(rbind, out)

  } else {
    if ((Sys.info()['sysname'] == "Windows")) {
      verbosecat("> SKAT processing multicore with snowparam on windows.\n")

      param  = BiocParallel::SnowParam(workers = n.cores, exportglobals = FALSE)
      wrap_fun = function(x, fx) {suppressPackageStartupMessages({
        library(dplyr); library(tibble)})
        fx(icr_ids[x], cpg_betas, df_study, model_str,
           cpg_mapping, m_value_transform, scaling = scaling, method = method, out_type = out_type)
      }
      out <- BiocParallel::bplapply(X = 1:length(icr_ids), FUN = wrap_fun, fx = skat_single_icr, BPPARAM = param)
      df_results = do.call(rbind, out)

    } else {
      verbosecat("> SKAT processing multicore with multicoreparam on mac/linux.\n")
    }
  }

  # Calculate adjusted p-value and q-value
  df_results$skat_adj_p_value <- p.adjust(p = df_results$skat_p_value, method = "fdr")
  df_results$skat_q_value <- qvalue::qvalue(p = df_results$skat_p_value, fdr.level = 0.05)$qvalues
  df_results <- df_results %>% dplyr::arrange(skat_adj_p_value)

  verbosecat(sprintf("> Filtering out ICRs with < %d cpg sites...\n", min_cpg))
  df_results <- df_results%>% filter(n_cpg >= min_cpg)

  return(df_results)
}



skat_single_icr <- function(icr_id, cpg_betas, df_study, model_str, cpg_mapping, m_value_transform, scaling, method, out_type) {
  # Get list of cpgs for a given ICR
  subset_cpg_ids <- cpg_mapping %>%
    dplyr::filter(.data$icr_id == .env$icr_id) %>% dplyr::pull(cpg_id)

  # Subset cpg beta matrix to only those contained within ICR
  tZ <- cpg_betas %>% tibble::rownames_to_column("cpg_id") %>% dplyr::filter(cpg_id %in% subset_cpg_ids) %>%
    tibble::column_to_rownames("cpg_id")
  # Tranform to m-values if requested
  if (m_value_transform) tZ = sesame::BetaValueToMValue(tZ)
  # Transform for input into null model
  #     rows = samples; columns = cpg sites
  Zs = Matrix::t(as.matrix(tZ))
  # Scaling (Center data)
  if (scaling) Zs <- scale(Zs, center = TRUE, scale = TRUE)

  # Calculate SKAT null model
  skat_null <- SKAT::SKAT_Null_Model(
    formula =  stats::as.formula(model_str), data = df_study,
    out_type = out_type, n.Resampling = 0, Adjustment = TRUE)

  # SKAT observed model
  skat_out<-SKAT::SKAT(Z = Zs, obj = skat_null, kernel = "linear", method = method)

  # Export data
  out <- data.frame(icr_id = icr_id, skat_p_value = skat_out$p.value, n_cpg = nrow(tZ))
  return(out)
}
