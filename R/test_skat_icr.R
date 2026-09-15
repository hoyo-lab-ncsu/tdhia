

#' Test ICR Methylation Associations with SKAT
#'
#' Fits a covariate-only null model and tests the joint association of each
#' ICR's CpG measurements with a continuous or dichotomous response using
#' a linear SKAT kernel.
#'
#' @param cpg_betas Numeric data frame of beta values with CpG IDs in rows
#'   and sample IDs in columns. The per-ICR helper requires data-frame input.
#' @param df_study Data frame with response and predictor columns. Row names
#'   must match the sample IDs in cpg_betas.
#' @param response Character string naming the response column.
#' @param predictors Character vector of covariate names for the null model.
#'   The current formula construction requires a nonempty predictor vector.
#' @param method Method passed to SKAT::SKAT(), such as "davies" or
#'   "optimal.adj" (the default SKAT-O combination of SKAT and burden tests).
#' @param out_type Response type passed to SKAT::SKAT_Null_Model():
#'   "C" for continuous (default) or "D" for dichotomous.
#' @param icr_ids Character vector of ICR IDs. NULL selects regions covered
#'   by the supplied CpGs in manifest_v1A2_design_scores.
#' @param min_cpg Minimum number of CpGs for retaining an ICR. Regions are
#'   tested first, then filtered before p-value adjustment.
#' @param db_flag Logical; save the initial environment to skat_icr_test.RData
#'   in the working directory. Defaults to FALSE.
#' @param m_value_transform Logical; convert beta values to M-values within
#'   each ICR before testing.
#' @param scaling Logical; center and scale each CpG across samples.
#' @param verbose Logical; print progress and the null-model formula.
#' @param n.cores Number of workers. One runs sequentially; larger values use
#'   SnowParam on Windows. The non-Windows multicore branch currently does
#'   not execute tests; use one core on those systems.
#'
#' @details
#' Samples with missing response or predictor values are removed, and the
#' remaining study rows are aligned to methylation columns by sample ID.
#' This function does not impute missing CpG measurements. The null model
#' uses Adjustment = TRUE and n.Resampling = 0.
#'
#' @return A data frame sorted by skat_adj_pvalue, containing icr_id,
#'   skat_raw_pvalue, n_cpg, skat_adj_pvalue (Benjamini-Hochberg adjustment),
#'   and skat_qvalue (from qvalue::qvalue()). Adjustments are calculated
#'   across regions remaining after the min_cpg filter.
#' @seealso [pc_regression_test()], [tdhia_stat_tests()]
#' @author Kate Everly, Bruce Corliss
#' @export
skat_icr_test <- function(cpg_betas, df_study, response, predictors,
                          method = "optimal.adj",
                          out_type="C", icr_ids = NULL,
                          min_cpg = 3, db_flag = FALSE, m_value_transform = T,
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
  df_study <- df_study %>% dplyr::select(dplyr::all_of(c(response, predictors))) %>% tidyr::drop_na()

  # Make the samples match and order them the same
  verbosecat("> Forcing sample id order to match between data and df_study.\n")
  shared_sample_ids <- intersect(colnames(cpg_betas),  rownames(df_study))
  cpg_betas <- cpg_betas[,shared_sample_ids]
  df_study <- df_study[rownames(df_study) %in% shared_sample_ids, ]

  # Order df_study the same sample order as cpg_beta
  df_study <- df_study[match(x = colnames(cpg_betas), table = rownames(df_study)),]

  # Dataframe to store results of test
  df_results <- data.frame(icr_id = icr_ids, skat_raw_pvalue = NA, n_cpg = NA)

  model_str = paste0(response , " ~ ", paste(predictors, collapse = " + "))
  verbosecat(sprintf("> Model: %s\n", model_str))

  
  # Calculate SKAT null model
  skat_null <- SKAT::SKAT_Null_Model(
    formula =  stats::as.formula(model_str), data = df_study,
    out_type = out_type, n.Resampling = 0, Adjustment = TRUE)
  
  # packaged skat into function for single ICR to make parallel conversion easier in future
  if (n.cores==1) { verbosecat("> SKAT processing on single core.\n")

    out = list()
    for (n in 1:length(icr_ids)) {

      out[[n]] <- skat_single_icr(
        icr_id = icr_ids[n], cpg_betas = cpg_betas, df_study = df_study, skat_null = skat_null,
        cpg_mapping = cpg_mapping, m_value_transform = m_value_transform, scaling = scaling, 
        method = method, out_type = out_type)
    }; df_results = do.call(rbind, out)

  } else {
    if ((Sys.info()['sysname'] == "Windows")) {
      verbosecat("> SKAT processing multicore with snowparam on windows.\n")

      param  = BiocParallel::SnowParam(workers = n.cores, exportglobals = FALSE)
      wrap_fun = function(x, fx) {suppressPackageStartupMessages({
        requireNamespace("dplyr"); requireNamespace("tibble")})
        fx(icr_ids[x], cpg_betas = cpg_betas, df_study = df_study, skat_null = skat_null,
           cpg_mapping = cpg_mapping, m_value_transform = m_value_transform, scaling = scaling, 
           method = method, out_type = out_type)
      }
      out <- BiocParallel::bplapply(X = 1:length(icr_ids), FUN = wrap_fun, fx = skat_single_icr, BPPARAM = param)
      df_results = do.call(rbind, out)

    } else {
      verbosecat("> SKAT processing multicore with multicoreparam on mac/linux.\n")
    }
  }

  # Remove ICRs with too few cpg sites
  verbosecat(sprintf("> Filtering out ICRs with < %d cpg sites...\n", min_cpg))
  df_results <- df_results %>% dplyr::filter(n_cpg >= min_cpg)
  
  # Calculate adjusted p-value and q-value
  df_results$skat_adj_pvalue <- p.adjust(p = df_results$skat_raw_pvalue, method = "fdr")
  df_results$skat_qvalue <- qvalue::qvalue(p = df_results$skat_raw_pvalue, fdr.level = 0.05)$qvalues
  df_results <- df_results %>% dplyr::arrange(skat_adj_pvalue)

  return(df_results)
}



skat_single_icr <- function(icr_id, cpg_betas, df_study, skat_null, cpg_mapping, m_value_transform, scaling, method, out_type) {
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

  # # Calculate SKAT null model
  # skat_null <- SKAT::SKAT_Null_Model(
  #   formula =  stats::as.formula(model_str), data = df_study,
  #   out_type = out_type, n.Resampling = 0, Adjustment = TRUE)

  # SKAT observed model
  skat_out<-SKAT::SKAT(Z = Zs, obj = skat_null, kernel = "linear", method = method)

  # Export data
  out <- data.frame(icr_id = icr_id, skat_raw_pvalue = skat_out$p.value, n_cpg = nrow(tZ))
  return(out)
}
