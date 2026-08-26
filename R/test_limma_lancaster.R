

#' Test CpG-wise differential methylation using limma
#'
#' @description
#' Tests associations between CpG methylation and study variables using
#' CpG-wise linear models fitted with \code{limma}.
#'
#' By default, input beta values are transformed to M-values before model
#' fitting because M-values generally have more appropriate variance properties
#' for linear-model inference. The original beta values are retained for
#' calculation of interpretable methylation effect sizes.
#'
#' The fitted model has the general form:
#'
#' \deqn{methylation \sim primary\ predictor + covariates}
#'
#' When \code{beadchip_correction = TRUE}, the variables \code{Bead},
#' \code{Col}, and \code{Row} are added to the model:
#'
#' \deqn{methylation \sim primary\ predictor + covariates + Bead + Col + Row}
#'
#' @details
#' The first entry in \code{predictors} is treated as the primary predictor of
#' interest. All remaining entries are included as adjustment covariates.
#'
#' The primary predictor must correspond to exactly one column in the model
#' matrix. This generally means it must be either:
#'
#' \itemize{
#'   \item a numeric variable, or
#'   \item a factor with exactly two levels.
#' }
#'
#' Predictors that expand to multiple model coefficients, such as factors with
#' more than two levels, are not currently supported unless an explicit contrast
#' is constructed outside this function.
#'
#' When the primary predictor is a two-level factor, the function calculates:
#'
#' \itemize{
#'   \item the mean beta value in the reference group,
#'   \item the mean beta value in the comparison group, and
#'   \item \code{delta_beta}, defined as comparison minus reference.
#' }
#'
#' A numeric primary predictor does not have natural reference and comparison
#' groups, so \code{delta_beta} is not calculated.
#'
#' When \code{m_value_transform = TRUE}, the \code{logFC} column returned by
#' \code{limma::topTable()} is on the M-value scale. It must not be interpreted
#' as a percentage methylation difference. For two-level factors,
#' \code{delta_beta} provides the corresponding effect size on the beta-value
#' scale.
#'
#' Samples are restricted to those present in both \code{df_study} and
#' \code{cpg_beta}. The methylation matrix is reordered to match the retained
#' metadata rows.
#'
#' CpGs are labeled \code{"Hyper"} or \code{"Hypo"} according to the sign of
#' \code{logFC} when the raw limma p-value is below
#' \code{pvalue_threshold}. Otherwise, they are labeled \code{"no"}.
#'
#' @param df_study A data frame containing one row per sample and columns for
#'   the sample identifier, primary predictor, covariates, and any requested
#'   beadchip adjustment variables. Categorical variables should be represented
#'   as factors with the intended reference level set before calling the
#'   function.
#'
#' @param predictors A character vector naming columns in \code{df_study} to
#'   include in the limma model. The first entry is the primary predictor tested;
#'   subsequent entries are adjustment covariates.
#'
#' @param cpg_beta A numeric matrix of methylation beta values, with CpG probes
#'   in rows and samples in columns. Row names must identify CpG probes. Column
#'   names must match sample identifiers in \code{df_study[[sample_name]]}.
#'
#' @param pvalue_threshold A numeric scalar between 0 and 1 giving the raw
#'   limma p-value threshold used to label CpGs as hypermethylated or
#'   hypomethylated and to select CpGs for annotation and plot labeling.
#'   Default is \code{0.0001}.
#'
#' @param db_flag Logical. If \code{TRUE}, saves the current function
#'   environment to \file{cpg_dml_test.RData} for debugging. Default is
#'   \code{FALSE}.
#'
#' @param sample_name A character scalar naming the sample-identifier column in
#'   \code{df_study}. Values in this column must match the column names of
#'   \code{cpg_beta}. Default is \code{"Patient.ID"}.
#'
#' @param beadchip_correction Logical. If \code{TRUE}, adds \code{Bead},
#'   \code{Col}, and \code{Row} to the model if they are not already included
#'   in \code{predictors}. These columns must exist in \code{df_study}. Default
#'   is \code{TRUE}.
#'
#' @param verbose Logical. If \code{TRUE}, prints progress and model information
#'   to the console. Default is \code{TRUE}.
#'
#' @param write_plots Logical. If \code{TRUE}, writes generated plots to
#'   \code{output_dir_path}. Default is \code{FALSE}.
#'
#' @param m_value_transform Logical. If \code{TRUE}, transforms beta values to
#'   M-values using \code{sesame::BetaValueToMValue()} before fitting the limma
#'   model. Beta values are still used to calculate group means and
#'   \code{delta_beta}. If \code{FALSE}, limma is fitted directly to beta
#'   values. Default is \code{TRUE}.
#'
#' @param correlation_check Logical. If \code{TRUE}, calculates and plots
#'   pairwise correlations among model predictors after converting
#'   non-numeric predictors to numeric factor codes. Default is \code{FALSE}.
#'
#' @param output_dir_path A character scalar giving the directory in which plot
#'   files are written when \code{write_plots = TRUE}. Default is the current
#'   working directory.
#'
#' @param label A character scalar used in output plot filenames. By default,
#'   the name of the primary predictor is used.
#'
#' @return A named list containing:
#'
#' \describe{
#'   \item{\code{df_dml}}{
#'     A data frame containing the complete
#'     \code{limma::topTable()} output for all tested CpGs, CpG and ICR
#'     annotations, differential-methylation labels, and genomic plotting
#'     coordinates. For a two-level primary predictor, the table also contains
#'     group-specific mean beta values and \code{delta_beta}.
#'   }
#'   \item{\code{chr_lens}}{
#'     A data frame containing chromosome plotting ranges and cumulative
#'     offsets used to construct the Manhattan plot.
#'   }
#'   \item{\code{cpgs_show}}{
#'     An annotated data frame restricted to CpGs with raw p-values below
#'     \code{pvalue_threshold}. When \code{delta_beta} is available, it also
#'     contains \code{beta_difference_percent}, which is
#'     \code{100 * delta_beta} formatted as an absolute percentage-point
#'     difference.
#'   }
#'   \item{\code{plots}}{
#'     A named list containing generated plot objects, including the p-value
#'     QQ plot, volcano plot, Manhattan plot, and optionally the predictor
#'     correlation plot.
#'   }
#' }
#'
#' @seealso
#' \code{\link[limma]{lmFit}},
#' \code{\link[limma]{eBayes}},
#' \code{\link[limma]{topTable}},
#' \code{\link[sesame]{BetaValueToMValue}},
#' \code{\link{icr_dmr_test}}
#'
#' @examples
#' \dontrun{
#' df_study$disease_state <- relevel(
#'   factor(df_study$disease_state),
#'   ref = "Control"
#' )
#'
#' dml_result <- cpg_dml_test(
#'   df_study = df_study,
#'   predictors = c("disease_state", "age", "sex"),
#'   cpg_beta = cpg_beta,
#'   sample_name = "Patient.ID",
#'   beadchip_correction = TRUE,
#'   m_value_transform = TRUE
#' )
#'
#' head(dml_result$df_dml)
#' }
#'
#' @export
cpg_dml_test <- function(df_study, predictors, cpg_beta,
                         pvalue_threshold = 0.0001, db_flag = F, sample_name = "Patient.ID",
                         beadchip_correction = T, verbose = T, write_plots = F,
                         m_value_transform = T,
                         correlation_check = F,
                         output_dir_path = getwd(), label = predictors[1]) {
  
  if (db_flag) {save(list = ls(all.names = TRUE), file = "cpg_dml_test.RData")}
  # load(file = "cpg_dml_test.RData")
  
  vsprintf <- \(x, ...) if (verbose) cat(sprintf(x, ...))
  plots = list()
  if (!sample_name %in% colnames(df_study)) {
    stop(sprintf(">> Error: sample_name column '%s' does not exist in df_study.", sample_name))
  }
  
  # Keep samples (patients) that exist in both study metadata and cpg_beta matrix
  shared_patients <- unique(df_study[[sample_name]][df_study[[sample_name]] %in% colnames(cpg_beta)])
  vsprintf(">> %i shared patients between study data and beta matrix, other culled.\n", length(shared_patients))
  if (length(shared_patients) == 0) {
    stop(">> Error: no shared patients between df_study and cpg_beta column names.")
  }
  df_study <- df_study[match(shared_patients, df_study[[sample_name]]), , drop = FALSE]
  cpg_beta <- cpg_beta[, shared_patients, drop = FALSE]
  
  
  # Transform cpg_beta
  # Use M-values for LIMMA inference when requested
  cpg_model_values <- if (m_value_transform) {
    sesame::BetaValueToMValue(cpg_beta)
  } else {
    cpg_beta
  }
  
  # Assemble predictor list
  if (beadchip_correction) {
    predictors <- unique(c(predictors, c("Bead", "Col", "Row")))
  }
  
  # Check that all required columns exist in study data.frame
  missed_columns <- setdiff(predictors, colnames(df_study))
  if (length(missed_columns)>0) {
    stop(sprintf(">> Error: these required columns do not exist in df_study: %s\n", 
             paste(missed_columns, collapse = ", ")))
  }
  if (anyNA(df_study[, predictors, drop = FALSE])) {
    stop(">> Error: df_study contains NA values in model predictor columns. Remove or impute them before cpg_dml_test().")
  }
  
  
  # Correlation plot of predictors                                    ##########
  #_____________________________________________________________________________
  if (correlation_check) {
    corr_data = df_study[, predictors] %>% 
      dplyr::mutate(dplyr::across(dplyr::everything(), ~ if (is.numeric(.x)) .x else as.numeric(as.factor(.x)))) %>% 
      as.data.frame()
    
    M <- stats::cor(corr_data, use = "pairwise.complete.obs")
    corr_result = corrplot::cor.mtest(corr_data, conf.level = 0.95)
    plots$predictor_correlations = corrplot::corrplot(
      M,p.mat = corr_result$p, insig = 'label_sig',
      sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9)
    
    if (write_plots) {
      grDevices::png(file.path(output_dir_path, paste0("corr_matrix_", label, ".png")),
          res = 300, width = 10, height = 10, units = "in")
      corrplot::corrplot(M,p.mat = corr_result$p, insig = 'label_sig',
                         sig.level = c(0.001, 0.01, 0.05),pch.cex = 0.9)
      grDevices::dev.off()
    }
  }
  
  # Get icr metadata
  ICR_CpG <- get_icr_metadata_table()
  
  
  # Create the model string to specify LIMMA analysis
  model_str = paste0("~ ", paste(predictors, collapse = " + "))
  vsprintf(">> Model string used for design matrix: %s\n", model_str)
  # Create design matrix based on model string
  model_formula <- stats::as.formula(model_str)
  design <- stats::model.matrix(model_formula, data = df_study)
  terms_obj <- stats::terms(model_formula)
  primary_term_index <- match(predictors[1], attr(terms_obj, "term.labels"))
  primary_coef_cols <- which(attr(design, "assign") == primary_term_index)
  if (length(primary_coef_cols) == 0) {
    stop(sprintf(">> Error: could not identify LIMMA coefficient for primary predictor '%s'.", predictors[1]))
  }
  if (length(primary_coef_cols) > 1) {
    stop(sprintf(
      ">> Error: primary predictor '%s' maps to multiple design columns (%s). Use a two-level factor, numeric predictor, or define an explicit contrast before calling cpg_dml_test().",
      predictors[1], paste(colnames(design)[primary_coef_cols], collapse = ", ")
    ))
  }
  coef_name <- colnames(design)[primary_coef_cols]
  non_estimable <- limma::nonEstimable(design)
  if (!is.null(non_estimable)) {
    stop(sprintf(">> Error: design matrix has non-estimable coefficients: %s", paste(non_estimable, collapse = ", ")))
  }
  vsprintf(">> LIMMA coefficient tested: %s\n", coef_name)

  # Determine reference and comparison groups
  #_____________________________________________________________________________
  primary_predictor <- predictors[1]
  primary_variable <- df_study[[primary_predictor]]
  if (all(primary_variable %in% 0:1)) primary_variable <- as.factor(primary_variable)
  beta_effects <- NULL
  
  if (is.factor(primary_variable) && nlevels(primary_variable) == 2L) {
    
    group_levels <- levels(primary_variable)
    reference_level <- group_levels[1]
    comparison_level <- group_levels[2]
    
    reference_samples <- df_study[[sample_name]][primary_variable == reference_level]
    
    comparison_samples <- df_study[[sample_name]][
      primary_variable == comparison_level
    ]
    
    mean_beta_reference <- rowMeans(
      cpg_beta[, reference_samples, drop = FALSE], na.rm = TRUE)
    
    mean_beta_comparison <- rowMeans(
      cpg_beta[, comparison_samples, drop = FALSE], na.rm = TRUE)
    
    beta_effects <- data.frame(
      CpG_Probe = rownames(cpg_beta),
      mean_beta_reference = mean_beta_reference,
      mean_beta_comparison = mean_beta_comparison,
      delta_beta = mean_beta_comparison - mean_beta_reference,
      stringsAsFactors = FALSE
    )
    
    names(beta_effects)[2:3] <- c(
      paste0("mean_beta_", make.names(reference_level)),
      paste0("mean_beta_", make.names(comparison_level))
    )
  }
  
  # Limma Analysis                                                  ############
  #_____________________________________________________________________________
  fit = limma::lmFit(cpg_model_values, design)
  fit = limma::eBayes(fit)
  df_dml = limma::topTable(fit, adjust.method="fdr", number=nrow(cpg_model_values), coef=coef_name)
  df_dml$CpG_Probe = rownames(df_dml)
  
  if (!is.null(beta_effects)) {
    df_dml <- df_dml %>%
      dplyr::left_join(beta_effects, by = "CpG_Probe")
  }
  
  
  df_dml = df_dml %>% dplyr::left_join(ICR_CpG, by = dplyr::join_by(CpG_Probe == CpG_id)) %>% 
    dplyr::arrange(adj.P.Val)
  df_dml$diffMeth = "no"
  df_dml$diffMeth[df_dml$P.Value < pvalue_threshold & df_dml$logFC > 0] = "Hyper"
  df_dml$diffMeth[df_dml$P.Value < pvalue_threshold & df_dml$logFC < 0] = "Hypo"
  
  # save(df_dml, file=paste0("EWAS/df_dml_EWAS_", label, ".Rdata"))
  # write.csv(df_dml, paste0("EWAS/df_dml_EWAS_", label, ".csv"))
  
  

  # QQ Plot Function
  #____________________________________________________________________________
  gg_qqplot = function(ps, ci = 0.95) {
    n  = length(ps)
    df = data.frame(
      observed = -log10(sort(ps)),
      expected = -log10(stats::ppoints(n)),
      clower   = -log10(stats::qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
      cupper   = -log10(stats::qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
    )
    log10Pe = expression(paste("Expected -log"[10], plain(P)))
    log10Po = expression(paste("Observed -log"[10], plain(P)))
    ggplot(df) +
      ggplot2::geom_ribbon(
        mapping = ggplot2::aes(x = expected, ymin = clower, ymax = cupper),
        alpha = 0.1
      ) +
      ggplot2::geom_point(ggplot2::aes(expected, observed), shape = 1, size = 3) +
      ggplot2::geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
      ggplot2::xlab(log10Pe) +
      ggplot2::ylab(log10Po)
  }
  
  inflation = function(ps) {
    chisq = stats::qchisq(1 - ps, 1)
    lambda = stats::median(chisq) / stats::qchisq(0.5, 1)
    lambda
  }
  
  
  
  # QQ and Volcano Plots
  #____________________________________________________________________________
  plots$qq_p_value = gg_qqplot(df_dml$P.Value) +
    ggplot2::theme_bw(base_size = 24) +
    ggplot2::annotate(geom = "text", x = -Inf, y = Inf, hjust = -0.15, vjust = 1 + 0.15 * 3,
                      label = sprintf("N; = %.2f", inflation(df_dml$P.Value)),
             size = 8) +
    ggplot2::theme(axis.ticks = element_line(linewidth = 0.5),
          panel.grid = element_blank())
  
  if (write_plots) {
    grDevices::png(file.path(output_dir_path, paste0("qqplot_", label, ".png")), res = 300, width = 10, height = 10, units = "in")
    print(plots$qq_p_value)
    grDevices::dev.off()
  }
  
  plots$volcano = ggplot2::ggplot(df_dml, ggplot2::aes(x=logFC, y = -log10(P.Value), col = diffMeth)) +
    ggplot2::geom_point() + ggplot2::theme_minimal() +
    ggrepel::geom_text_repel(data = df_dml[df_dml$P.Value < pvalue_threshold,], ggplot2::aes(
      label=paste0(CpG_Probe, " (", ICR_id , ")")), size=3  )
  
  if (write_plots) {
    grDevices::png(file.path(output_dir_path, paste0("VolcanoPlot_", label, ".png")), res = 300, width = 12, height = 8, units = "in")
    print(plots$volcano)
    grDevices::dev.off()
  }
  
  # Manhattan Plots
  # ____________________________________________________________________________
  
  chr_lens <- df_dml %>% dplyr::group_by(chr_num) %>% 
    dplyr::summarise(chr_len = as.double(max(CpG_start)))
  chr_lens <- chr_lens %>% dplyr::mutate(offset = dplyr::lag(cumsum(chr_len), 
                                                             default = 0))
  
  df_dml <- dplyr::inner_join(df_dml, chr_lens, by = "chr_num") %>%
    dplyr::mutate(bp_cum = CpG_start + offset)
  axis_df <- df_dml %>% dplyr::group_by(chr_num) %>% 
    dplyr::summarise(center = mean(bp_cum))
  # sig <- sig
  
  plots$manhattan = ggplot2::ggplot(df_dml, ggplot2::aes(x = bp_cum, y = -log10(P.Value), color = diffMeth)) +
    ggplot2::geom_point(alpha = 0.75) + ggplot2::geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed") +
    ggplot2::scale_x_continuous(breaks = axis_df$center, labels = axis_df$chr_num)  +
    ggplot2::scale_color_manual(values = c(Hyper = "#F8766D", Hypo = "#35C96E", no = "grey70"))+
    ggplot2::theme_minimal() +
    ggplot2::theme( 
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      axis.title.y = ggtext::element_markdown(),
      axis.text.x = ggplot2::element_text(angle = 60, size = 8, vjust = 0.5)
    )+
    labs(y = "-log10(p)", x = "Chromosome")+ 
    ggrepel::geom_text_repel(
      data = df_dml[which(df_dml$P.Value<pvalue_threshold),],
      ggplot2::aes(x=bp_cum, y = -log10(P.Value), label=paste0(CpG_Probe, " (", ICR_id , ")")),
      box.padding = 0.5, point.padding = 0.3, max.overlaps = 10, size = 3.5, nudge_x = 0.05, 
      nudge_y = 0.05)
  
  if (write_plots) {
    grDevices::png(file.path(output_dir_path, paste0("Manhattan_", label, ".png")), res = 300, width = 15, height = 5, units = "in")
    print(plots$manhattan)
    grDevices::dev.off()
  }
  

  # Gene annotation for significant cpgs                    ####################
  #_____________________________________________________________________________
  
  cpgs <- df_dml %>% dplyr::filter(P.Value < pvalue_threshold) %>% dplyr::arrange(P.Value)
  
  # Get gene info
  genes_gr <- get_gene_info_biomart()
  
  
  promoters_gr <- GenomicRanges::promoters(genes_gr, upstream = 1500, downstream = 500)
  
  cpg_gr <- GenomicRanges::GRanges(seqnames = cpgs$CpG_chr,
                    ranges = IRanges::IRanges(start = cpgs$CpG_start, end = cpgs$CpG_stop),
                    CpG_Probe = cpgs$CpG_Probe)
  
  overlaps <- GenomicRanges::findOverlaps(cpg_gr, genes_gr)
  
  genes_overlap <- data.frame(CpG_Probe = cpg_gr$CpG_Probe[S4Vectors::queryHits(overlaps)],
                              Gene = genes_gr$gene_name[S4Vectors::subjectHits(overlaps)])
  
  cpg_gr_extended <- GenomicRanges::resize(cpg_gr, width = 20001, fix = "center")
  
  nearby <- GenomicRanges::findOverlaps(cpg_gr_extended, genes_gr)
  genes_nearby <- data.frame(CpG_Probe = cpg_gr$CpG_Probe[S4Vectors::queryHits(nearby)],
                             Gene = genes_gr$gene_name[S4Vectors::subjectHits(nearby)])
  
  all_annotations <- dplyr::bind_rows(genes_overlap, genes_nearby) %>% 
    dplyr::distinct(CpG_Probe, Gene) %>% dplyr::group_by(CpG_Probe) %>% 
    dplyr::summarise(Gene = paste(unique(Gene), collapse = ", "), .groups = "drop")
  
  cpgs_annotated <- merge(cpgs, all_annotations, by = "CpG_Probe", all.x = TRUE)
  
  cpgs_show <- cpgs_annotated 
  
  
  
  if ("delta_beta" %in% names(cpgs_show)) {
    cpgs_show <- cpgs_show %>%
      dplyr::mutate(
        beta_difference_percent = sprintf(
          "%+.2f%%",
          delta_beta * 100
        )
      ) %>%
      dplyr::relocate(
        delta_beta,
        beta_difference_percent,
        .after = logFC
      )
  }

  
  # Export figures and data
  out <- list(df_dml = df_dml, chr_lens = chr_lens, cpgs_show = cpgs_show, 
              plots = plots)
  
  return(out)
}




#' Aggregate CpG-level methylation results by imprinting control region
#'
#' @description
#' Performs an imprinting control region-level analysis by aggregating
#' CpG-level differential methylation results produced by
#' \code{cpg_dml_test()}.
#'
#' CpGs are grouped by imprinting control region (ICR). Within each ICR,
#' CpG-level p-values are combined using the Lancaster method, with the absolute
#' CpG-level \code{logFC} values used as weights.
#'
#' @details
#' The function summarizes each ICR using:
#'
#' \itemize{
#'   \item the number of represented CpGs,
#'   \item the mean CpG-level \code{logFC},
#'   \item the mean \code{delta_beta}, when available,
#'   \item the number of CpGs with positive and negative coefficients,
#'   \item the percentage of CpGs with raw p-values below
#'     \code{pval_threshold},
#'   \item a Lancaster combined p-value, and
#'   \item a Benjamini-Hochberg FDR-adjusted combined p-value.
#' }
#'
#' The regional direction is labeled \code{"Hyper"} when more CpGs have
#' positive than negative \code{logFC} values. Otherwise, the region is labeled
#' \code{"Hypo"}. Ties are therefore labeled \code{"Hypo"} by the current
#' implementation.
#'
#' When CpG-level models were fitted to M-values, \code{mean_logFC} is also on
#' the M-value scale and must not be interpreted as a percentage methylation
#' difference. When available, \code{mean_delta_beta} provides the average
#' CpG-level difference on the beta-value scale.
#'
#' ICRs are annotated using overlapping, promoter-associated, and nearby
#' protein-coding genes. Nearby genes are identified using an ICR interval
#' extended by 10,000 base pairs in each direction.
#'
#' The function also creates a regional Manhattan plot using the chromosome
#' offsets supplied in \code{chr_lens}.
#'
#' @param df_dml A data frame containing CpG-level differential methylation
#'   results, normally the \code{df_dml} component returned by
#'   \code{cpg_dml_test()}.
#'
#'   The data frame must contain at least:
#'
#'   \itemize{
#'     \item \code{ICR_id}: imprinting control region identifier;
#'     \item \code{ICR_chr}: ICR chromosome;
#'     \item \code{ICR_start}: ICR start coordinate;
#'     \item \code{ICR_stop}: ICR end coordinate;
#'     \item \code{logFC}: CpG-level limma coefficient;
#'     \item \code{P.Value}: raw CpG-level p-value; and
#'     \item \code{adj.P.Val}: FDR-adjusted CpG-level p-value.
#'   }
#'
#'   An optional \code{delta_beta} column may contain beta-value differences
#'   for a two-level primary predictor. If it is absent, the function adds it
#'   internally with missing values.
#'
#' @param chr_lens A data frame of chromosome plotting ranges, normally the
#'   \code{chr_lens} component returned by \code{cpg_dml_test()}. It must
#'   contain:
#'
#'   \itemize{
#'     \item \code{chr_num}: numeric chromosome identifier;
#'     \item \code{chr_len}: plotting extent for the chromosome; and
#'     \item \code{offset}: cumulative genomic-coordinate offset.
#'   }
#'
#' @param pval_threshold A numeric scalar between 0 and 1. Within each ICR, the
#'   percentage of significant CpGs is calculated as the percentage with raw
#'   \code{P.Value < pval_threshold}. Default is \code{0.05}.
#'
#' @param fdr_sig_threshold A numeric scalar between 0 and 1 giving the
#'   ICR-level FDR threshold used to label regions in the Manhattan plot.
#'   Default is \code{0.0001}.
#'
#' @param verbose Logical. Reserved for progress reporting. Default is
#'   \code{TRUE}.
#'
#' @param db_flag Logical. If \code{TRUE}, saves the current function
#'   environment to \file{icr_dmr_test.RData} for debugging. Default is
#'   \code{FALSE}.
#'
#' @return A named list containing:
#'
#' \describe{
#'   \item{\code{ICR_summary}}{
#'     One row per ICR, containing regional coordinates, CpG count,
#'     \code{mean_logFC}, \code{mean_delta_beta}, combined p-value, FDR,
#'     direction, counts of positive and negative CpG coefficients, and the
#'     percentage of CpGs passing \code{pval_threshold}.
#'   }
#'   \item{\code{ICR_collapsed}}{
#'     A gene-annotated table of prioritized ICRs after application of the
#'     regional filtering criteria defined in the function.
#'   }
#'   \item{\code{plots}}{
#'     A named list containing the regional Manhattan plot in
#'     \code{manhplot_icr}.
#'   }
#' }
#'
#' @seealso
#' \code{\link{cpg_dml_test}},
#' \code{\link[aggregation]{lancaster}}
#'
#' @examples
#' \dontrun{
#' dmr_result <- icr_dmr_test(
#'   df_dml = dml_result$df_dml,
#'   chr_lens = dml_result$chr_lens,
#'   pval_threshold = 0.05,
#'   fdr_sig_threshold = 0.0001
#' )
#'
#' head(dmr_result$ICR_summary)
#' }
#'
#' @export
icr_dmr_test <- function(df_dml, chr_lens, pval_threshold = 0.05,
                         fdr_sig_threshold = 0.0001,
                         verbose = T, db_flag = F) {
  if (db_flag) {save(list = ls(all.names = TRUE), file = "icr_dmr_test.RData")}
  # load(file = "icr_dmr_test.RData")
  

  # Regional (DMR-like) Analysis via ICR Aggregation
  # ____________________________________________________________________________


  
  data_signif <- df_dml %>% dplyr::filter(adj.P.Val < pval_threshold)
  
  # genes <- biomaRt::getBM(attributes = c("chromosome_name", "start_position", "end_position", "hgnc_symbol", "gene_biotype", "strand"),
  #                filters = "biotype", values = "protein_coding", mart = mart)
  # 
  # genes_gr <- GenomicRanges::GRanges(seqnames = paste0("chr", genes$chromosome_name),
  #                     ranges = IRanges::IRanges(start = genes$start_position, end = genes$end_position),
  #                     gene_name = genes$hgnc_symbol, gene_biotype = genes$gene_biotype, strand = genes$strand)
  
  genes_gr <- get_gene_info_biomart()
  
  promoters_gr <- GenomicRanges::promoters(genes_gr, upstream = 1500, downstream = 500)
  
  
  if (!"delta_beta" %in% names(df_dml)) {
    df_dml$delta_beta <- NA_real_
  }
  
  
  ICR_summary <- df_dml %>% 
    dplyr::group_by(ICR_id ) %>% 
    dplyr::summarise(ICR_n = dplyr::n(), 
              ICR_chr = unique(ICR_chr), 
              ICR_start = min(ICR_start), 
              ICR_stop = max(ICR_stop), 
              mean_logFC = mean(logFC, na.rm = TRUE), 
              mean_delta_beta = if (all(is.na(delta_beta))) {
                NA_real_
              } else {
                mean(delta_beta, na.rm = TRUE)
              },
              pval_combined = if (dplyr::n() > 1) {
                aggregation::lancaster(pvalues = P.Value, weights = abs(logFC))
              } else { 
                P.Value[1] 
              },
              direction = ifelse(sum(logFC > 0) > sum(logFC < 0), "Hyper", "Hypo"),
              n_hyper = sum(logFC > 0), n_hypo = sum(logFC < 0),
              perc_signif = (sum(P.Value < pval_threshold) / dplyr::n()) * 100) %>%
    dplyr::ungroup() %>% 
    dplyr::mutate(FDR = p.adjust(pval_combined, method = "fdr")) %>% 
    dplyr::arrange(pval_combined)
  
  ICR_gr <- GenomicRanges::GRanges(seqnames = ICR_summary$ICR_chr,
                    ranges = IRanges::IRanges(start = ICR_summary$ICR_start, end = ICR_summary$ICR_stop),
                    ICR_id = ICR_summary$ICR_id)
  
  overlaps <- IRanges::findOverlaps(ICR_gr, genes_gr)
  
  genes_overlap <- data.frame(ICR_id = ICR_gr$ICR_id[S4Vectors::queryHits(overlaps)],
                              Gene = genes_gr$gene_name[S4Vectors::subjectHits(overlaps)], Region = "Overlap")
  
  ICR_gr_extended <- GenomicRanges::resize(ICR_gr, GenomicRanges::width(ICR_gr) + 20000, fix = "center")
  
  proximal <- IRanges::findOverlaps(ICR_gr_extended, genes_gr)
  
  genes_nearby <- data.frame(ICR_id = ICR_gr$ICR_id[S4Vectors::queryHits(proximal)],
                             Gene = genes_gr$gene_name[S4Vectors::subjectHits(proximal)], Region = "Nearby")
  
  promoter_hits <- IRanges::findOverlaps(ICR_gr, promoters_gr)
  
  genes_promoter <- data.frame(ICR_id = ICR_gr$ICR_id[S4Vectors::queryHits(promoter_hits)],
                               Gene = promoters_gr$gene_name[S4Vectors::subjectHits(promoter_hits)], Region = "Promoter")
  
  all_annotations <- dplyr::bind_rows(genes_overlap, genes_nearby, genes_promoter) %>% 
    dplyr::distinct(ICR_id, Gene, .keep_all = TRUE)
  ICR_annotated <- merge(ICR_summary, all_annotations, by = "ICR_id", all.x = TRUE)
  
  ICR_annotated$Gene[ICR_annotated$Gene == ""] <- NA

  
  ICR_annotated_plot <- ICR_annotated %>%
    dplyr::group_by(ICR_id, ICR_n, ICR_chr, ICR_start, ICR_stop, mean_logFC,
                    mean_delta_beta, direction, n_hyper, n_hypo,  perc_signif, 
                    FDR ) %>%
    dplyr::summarise( Gene = paste(stats::na.omit(unique(Gene)), collapse = ", "),
      .groups = "drop") %>%
    dplyr::arrange(FDR) %>%
    dplyr::mutate(  perc_signif = paste0(round(perc_signif, 1), "%"),
      beta_difference_percent = dplyr::if_else(  is.na(mean_delta_beta),
        NA_character_, sprintf("%+.2f%%", mean_delta_beta * 100)  )) %>%
    dplyr::relocate(  mean_delta_beta,  beta_difference_percent,
      .after = mean_logFC )
  
  ICR_annotated_plot$ICR_chr[ICR_annotated_plot$ICR_chr == "chrX"] <- "chr23"
  ICR_annotated_plot$chr_num <- as.numeric(stringr::str_remove(ICR_annotated_plot$ICR_chr, "chr"))
  ICR_annotated_plot$start <- ICR_annotated_plot$ICR_start
  
  chr_lens <- chr_lens %>% dplyr::rename_with(~ "bp_add", .cols = "offset")
  ICR_annotated_plot <- ICR_annotated_plot %>% 
    dplyr::inner_join(chr_lens, by = "chr_num") %>% 
    dplyr::mutate(bp_cum = start + bp_add)
  
  axis_df <- ICR_annotated_plot %>% 
    dplyr::group_by(chr_num) %>% 
    dplyr::summarise(center = mean(bp_cum))
  
  ylim <- max(-log10(ICR_annotated_plot$FDR), na.rm = TRUE) + 1
  
  manhplot_icr <- ggplot(ICR_annotated_plot, aes(x = bp_cum, y = -log10(FDR), color = direction)) +
    ggplot2::geom_point(alpha = 0.75) +
    ggplot2::geom_hline(yintercept = -log10(0.0001), linetype = "dashed") +
    ggplot2::scale_x_continuous(breaks = axis_df$center, labels = axis_df$chr_num) +
    ggplot2::theme_minimal() +
    ggplot2:: theme( 
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      axis.title.y = ggtext::element_markdown(),
      axis.text.x = ggplot2::element_text(angle = 60, size = 8, vjust = 0.5)
    )+ 
    ggplot2::labs(x = NULL, y = "-log10(FDR)")+ 
    ggrepel::geom_text_repel(
      data = ICR_annotated_plot[which(ICR_annotated_plot$FDR<fdr_sig_threshold),],aes(x=bp_cum, y = -log10(FDR), label=ICR_id),
      box.padding = 0.5,
      point.padding = 0.3,
      max.overlaps = 10,
      size = 3.5,
      nudge_x = 0.05, 
      nudge_y = 0.05)
  
  
  # png(paste0(output_dir_path, "/Manhattan_ICR_", label, ".png"), res = 300, width = 15, height = 5, units = "in")
  # print(manhplot_icr)
  # dev.off()
  
  v_icr <- ggplot2::ggplot(ICR_annotated_plot, ggplot2::aes(x = mean_logFC, y = -log10(FDR), color = direction, label = ICR_id)) +
    ggplot2::geom_point() +
    ggrepel::geom_text_repel(data = ICR_annotated_plot %>% dplyr::filter(FDR < 0.0001), ggplot2::aes(label = ICR_id)) +
    ggplot2::theme_minimal()
  
  # png(paste0(output_dir_path, "/Volcano_ICR_", label, ".png"), res = 300, width = 12, height = 8, units = "in")
  # print(v_icr)
  # dev.off()
  
  
  ICR_prioritized <- ICR_annotated %>%
    dplyr::filter(
      FDR < 0.05,
      abs(mean_logFC) > 0.015,
      perc_signif > 40
    )
  
  ICR_prioritized$Gene[ICR_prioritized$Gene == ""] <- NA
  
  ICR_collapsed <- ICR_prioritized %>%
    dplyr::group_by(ICR_id, ICR_n, ICR_chr, ICR_start, ICR_stop, mean_logFC,
                    mean_delta_beta, direction, n_hyper, n_hypo,  perc_signif, 
                    FDR ) %>%
    dplyr::summarise( Gene = paste(stats::na.omit(unique(Gene)), collapse = ", "),
                      .groups = "drop") %>%
    dplyr::arrange(FDR) %>%
    dplyr::mutate(  perc_signif = paste0(round(perc_signif, 1), "%"),
                    beta_difference_percent = dplyr::if_else(  is.na(mean_delta_beta),
                                                               NA_character_, sprintf("%+.2f%%", mean_delta_beta * 100)  )) %>%
    dplyr::relocate(  mean_delta_beta,  beta_difference_percent,
                      .after = mean_logFC )
  
  
  out <- list(ICR_summary = ICR_summary, ICR_collapsed = ICR_collapsed,
              plots= list(manhplot_icr = manhplot_icr))
  # save(ICR_collapsed, file = paste0(output_dir_path, "/res_EWAS_", label, "_ICR.Rdata"))
  
  return(out)
  
  
}




#' Retrieve and cache human protein-coding gene annotations
#'
#' @description
#' Retrieves human protein-coding gene annotations from Ensembl BioMart and
#' returns them as a \code{GenomicRanges::GRanges} object.
#'
#' @details
#' On the first call, the function queries the
#' \code{ENSEMBL_MART_ENSEMBL} mart using the
#' \code{hsapiens_gene_ensembl} dataset. It retrieves chromosome, start
#' position, end position, HGNC symbol, and gene biotype for protein-coding
#' genes.
#'
#' The resulting \code{GRanges} object is saved as
#' \file{biomart_genelist.rds} in the package-specific user data directory
#' returned by:
#'
#' \preformatted{
#' tools::R_user_dir(utils::packageName(), which = "data")
#' }
#'
#' Subsequent calls load the cached object instead of querying BioMart again.
#' This reduces execution time and avoids repeated requests that may trigger
#' connection timeouts or server-side request limits.
#'
#' Chromosome names are prefixed with \code{"chr"}. Gene symbols are stored in
#' the \code{gene_name} metadata column, and gene biotypes are stored in the
#' \code{gene_biotype} metadata column.
#'
#' Because the BioMart result is cached without an expiration date, the returned
#' annotations may not reflect later Ensembl updates unless the cached file is
#' manually removed.
#'
#' @return A \code{GenomicRanges::GRanges} object containing human
#'   protein-coding gene intervals. Metadata columns include
#'   \code{gene_name} and \code{gene_biotype}.
#'
#' @seealso
#' \code{\link[biomaRt]{useMart}},
#' \code{\link[biomaRt]{getBM}},
#' \code{\link[tools]{R_user_dir}}
#'
#' @examples
#' \dontrun{
#' genes_gr <- get_gene_info_biomart()
#' genes_gr
#' S4Vectors::mcols(genes_gr)
#' }
#'
#' @export
get_gene_info_biomart <- function() {
  
  internal_package_path <- tools::R_user_dir(utils::packageName(), which = "data")
  biomart_path <- paste0(internal_package_path, "/biomart_genelist.rds")
  dir.create(internal_package_path, showWarnings = F, recursive = T)
  if ( !file.exists(biomart_path) ) {
    mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                             dataset = "hsapiens_gene_ensembl")
    #host = "https://useast.ensembl.org")
    genes <- biomaRt::getBM(attributes = c("chromosome_name", "start_position", "end_position", "hgnc_symbol", "gene_biotype"),
                            filters = "biotype", values = "protein_coding", mart = mart)
    
    genes_gr <- GenomicRanges::GRanges(seqnames = paste0("chr", genes$chromosome_name),
                                       ranges = IRanges::IRanges(start = genes$start_position, end = genes$end_position),
                                       gene_name = genes$hgnc_symbol, gene_biotype = genes$gene_biotype, strand = "*")
    saveRDS(genes_gr, file = biomart_path)
  } else {
    genes_gr <- readRDS(file = biomart_path)
  }
  
  
  
  return(genes_gr)
}


#' Retrieve CpG-to-ICR mapping metadata
#'
#' @description
#' Returns the package's CpG-to-imprinting-control-region mapping table and adds
#' a numeric chromosome identifier for plotting and chromosome-based
#' operations.
#'
#' @details
#' The function derives its output from
#' \code{tdhia::mapping_cpg_icr_ids}. The \code{"chr"} prefix is removed from
#' \code{ICR_chr} to create \code{chr_num}. Chromosome X is encoded as 23 and
#' chromosome Y is encoded as 24 before conversion to numeric.
#'
#' The original chromosome fields are retained unchanged.
#'
#' @return A data frame containing CpG probe and imprinting control region
#'   annotations. Expected columns include:
#'
#' \describe{
#'   \item{\code{CpG_chr}}{Chromosome containing the CpG probe.}
#'   \item{\code{CpG_start}}{CpG start coordinate.}
#'   \item{\code{CpG_stop}}{CpG end coordinate.}
#'   \item{\code{ICR_chr}}{Chromosome containing the associated ICR.}
#'   \item{\code{ICR_start}}{ICR start coordinate.}
#'   \item{\code{ICR_stop}}{ICR end coordinate.}
#'   \item{\code{CpG_Probe}}{CpG probe identifier, when present in the source
#'     table.}
#'   \item{\code{CpG_id}}{CpG identifier used to join the metadata to
#'     CpG-level limma results.}
#'   \item{\code{ICR_id}}{Imprinting control region identifier.}
#'   \item{\code{chr_num}}{Numeric ICR chromosome identifier, with X encoded
#'     as 23 and Y encoded as 24.}
#' }
#'
#' Additional columns present in \code{mapping_cpg_icr_ids} are retained.
#'
#' @examples
#' icr_metadata <- get_icr_metadata_table()
#' head(icr_metadata)
#'
#' @export
get_icr_metadata_table <- function() {
  ICR_CpG <- tdhia::mapping_cpg_icr_ids %>%
    dplyr::mutate(chr_num = stringr::str_replace(ICR_chr, "^chr", ""))
  ICR_CpG$chr_num[ICR_CpG$chr_num == "X"] <- 23
  ICR_CpG$chr_num[ICR_CpG$chr_num == "Y"] <- 24
  ICR_CpG$chr_num = as.numeric(ICR_CpG$chr_num)
  
  return(ICR_CpG)
}
