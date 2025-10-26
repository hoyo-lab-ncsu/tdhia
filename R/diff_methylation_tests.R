


#' cpg_dmr_test
#'
#' @description analyze for changes in beta value between a control group and an
#' experiment group.
#'
#' @param metadata a matrix of beta values, cpg/ ICR site (rows) by patients (columns)
#' @param betas boolean vector of which columns in beta matrix are control group
#' @export
#' @author author
cpg_dmr_test <- function(metadata, betas, analysis_type = "casecontrol",  sig = 0.0001, db_flag = T) {
  if (db_flag) {save(list = ls(all.names = TRUE), file = "cpg_dmr_test.RData")}
  # load(file = "cpg_dmr_test.RData")
  
  # -----------------------------
  # Annotation
  # -----------------------------
  
  ICR_CpG <- tdhia::mapping_cpg_icr_ids %>%
    dplyr::mutate(chr_num = stringr::str_replace(ICR_chr, "^chr", ""))
  ICR_CpG$chr_num[ICR_CpG$chr_num == "X"] <- 23
  ICR_CpG$chr_num[ICR_CpG$chr_num == "Y"] <- 24
  ICR_CpG$chr_num = as.numeric(ICR_CpG$chr_num)
  
  
  # -----------------------------
  # Design Matrix
  # -----------------------------
  
  rownames(metadata) <- metadata$Patient.ID
  metadata <- metadata[colnames(betas),]
  
  if (analysis_type == "casecontrol") {
    design <- stats::model.matrix(~ Group + Beadchip + Col + Row + Race, data = metadata)
    coef_name <- "GroupCase"
  } else {
    design <- stats::model.matrix(~ PCL.Total + Beadchip + Col + Row + Race, data = metadata)
    coef_name <- "PCL.Total"
  }
  design <- design[, colSums(design) != 1] 
  
  
  # -----------------------------
  # Limma Analysis
  # -----------------------------
  
  fit = limma::lmFit(betas, design)
  fit = limma::eBayes(fit)
  res = limma::topTable(fit, adjust.method="fdr", number=nrow(betas), coef=coef_name)
  res$CpG_Probe = rownames(res)
  res = res %>% dplyr::left_join(ICR_CpG, by=dplyr::join_by(CpG_Probe == CpG_id)) %>% 
    dplyr::arrange(adj.P.Val)
  res$diffMeth = "no"
  res$diffMeth[res$P.Value < 0.0001 & res$logFC > 0] = "Hyper"
  res$diffMeth[res$P.Value < 0.0001 & res$logFC < 0] = "Hypo"
  
  # save(res, file=paste0("EWAS/res_EWAS_", label, ".Rdata"))
  # write.csv(res, paste0("EWAS/res_EWAS_", label, ".csv"))
  
  
  # -----------------------------
  # QQ Plot Function
  # -----------------------------
  gg_qqplot = function(ps, ci = 0.95) {
    n  = length(ps)
    df = data.frame(
      observed = -log10(sort(ps)),
      expected = -log10(ppoints(n)),
      clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
      cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
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
    chisq = qchisq(1 - ps, 1)
    lambda = median(chisq) / qchisq(0.5, 1)
    lambda
  }
  
  
  
  # -----------------------------
  # QQ and Volcano Plots
  # -----------------------------
  p = gg_qqplot(res$P.Value) +
    ggplot2::theme_bw(base_size = 24) +
    ggplot2::annotate(geom = "text",
             x = -Inf,
             y = Inf,
             hjust = -0.15,
             vjust = 1 + 0.15 * 3,
             label = sprintf("N; = %.2f", inflation(res$P.Value)),
             size = 8) +
    ggplot2::theme(axis.ticks = element_line(linewidth = 0.5),
          panel.grid = element_blank())
  
  # png(paste0(output_dir_path, "qqplot_", label, ".png"), res = 300, width = 10, height = 10, units = "in")
  # print(p)
  # dev.off()
  # 
  res$diffMeth[res$logFC > 0] = "Hyper"
  res$diffMeth[res$logFC < 0] = "Hypo"
  
  volcano = ggplot2::ggplot(res, ggplot2::aes(x=logFC, y=-log10(P.Value), col=diffMeth)) +
    ggplot2::geom_point() + ggplot2::theme_minimal() +
    ggrepel::geom_text_repel(data = res[res$P.Value<0.0001,], ggplot2::aes(label=paste0(CpG_Probe, " (", ICR_id , ")")), size=3)
  
  # png(paste0(output_dir_path, "/VolcanoPlot_", label, ".png"), res = 300, width = 12, height = 8, units = "in")
  # print(volcano)
  # dev.off()
  
  
  # -----------------------------
  # Manhattan Plots
  # -----------------------------
  
  chr_lens <- res %>% dplyr::group_by(chr_num) %>% dplyr::summarise(chr_len = as.double(max(CpG_start)))
  chr_lens <- chr_lens %>% dplyr::mutate(offset = dplyr::lag(cumsum(chr_len), default = 0))
  res <- dplyr::inner_join(res, chr_lens, by = "chr_num") %>% dplyr::mutate(bp_cum = CpG_start + offset)
  axis_df <- res %>% dplyr::group_by(chr_num) %>% dplyr::summarise(center = mean(bp_cum))
  # sig <- 0.0001
  
  manhattan = ggplot(res, aes(x = bp_cum, y = -log10(P.Value), color = diffMeth)) +
    geom_point(alpha = 0.75) +
    geom_hline(yintercept = -log10(sig), linetype = "dashed") +
    scale_x_continuous(breaks = axis_df$center, labels = axis_df$chr_num)  +
    scale_color_manual(values = c("#F8766D", "#35C96E"))+
    theme_minimal() +
    theme( 
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = ggtext::element_markdown(),
      axis.text.x = ggplot2::element_text(angle = 60, size = 8, vjust = 0.5)
    )+
    labs(y = "-log10(p)", x = "Chromosome")+ 
    ggrepel::geom_text_repel(
      data = res[which(res$P.Value<sig),],aes(x=bp_cum, y = -log10(P.Value), label=paste0(CpG_Probe, " (", ICR_id , ")")),
      box.padding = 0.5,
      point.padding = 0.3,
      max.overlaps = 10,
      size = 3.5,
      nudge_x = 0.05, 
      nudge_y = 0.05)
  
  
  # png(paste0(output_dir_path, "/Manhattan_", label, ".png"), res = 300, width = 15, height = 5, units = "in")
  # print(manhattan)
  # dev.off()
  
  
  
  
  # -----------------------------
  # Gene Annotation for Significant CpGs
  # -----------------------------
  
  cpgs <- res %>% dplyr::filter(P.Value < 0.0001) %>% dplyr::arrange(P.Value)
  
  
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
  genes_nearby <- data.frame(CpG_Probe = cpg_gr$CpG_Probe[S4Vectors::queryHits(nearby)], Gene = genes_gr$gene_name[S4Vectors::subjectHits(nearby)])
  
  all_annotations <- dplyr::bind_rows(genes_overlap, genes_nearby) %>% 
    dplyr::distinct(CpG_Probe, Gene) %>% dplyr::group_by(CpG_Probe) %>% 
    dplyr::summarise(Gene = paste(unique(Gene), collapse = ", "), .groups = "drop")
  
  cpgs_annotated <- merge(cpgs, all_annotations, by = "CpG_Probe", all.x = TRUE)
  
  cpgs_show <- cpgs_annotated %>% 
    dplyr::mutate(perc_diff = ifelse(logFC >= 0, paste0("+", round(logFC * 100, 2), "%"), paste0("-", round(abs(logFC * 100), 2), "%"))) %>% 
    dplyr::relocate(perc_diff, .after = logFC)
  
  
  # Export figrues and data
  out <- list(res = res, chr_lens = chr_lens, cpgs_show = cpgs_show, 
              plots = list(manhattan = manhattan, p = p, volcano = volcano))
  
  return(out)
}




#' icr_dmr_test
#'
#' @description analyze for changes in beta value between a control group and an
#' experiment group.
#'
#' @param res a matrix of beta values, cpg/ ICR site (rows) by patients (columns)
#' @param chr_lens boolean vector of which columns in beta matrix are control group
#' @param pval_threshold max value for signficiant p-value from test (0-1 singleton).
#' @export
#' @author author
icr_dmr_test <- function(res, chr_lens, pval_threshold = 0.05, sig = 0.0001) {
  
  
  
  
  # -----------------------------
  # Regional (DMR-like) Analysis via ICR Aggregation
  # -----------------------------
  # pval_threshold <- 0.05
  res_pre <- res
  
  data_signif <- res_pre %>% dplyr::filter(adj.P.Val < pval_threshold)
  
  # genes <- biomaRt::getBM(attributes = c("chromosome_name", "start_position", "end_position", "hgnc_symbol", "gene_biotype", "strand"),
  #                filters = "biotype", values = "protein_coding", mart = mart)
  # 
  # genes_gr <- GenomicRanges::GRanges(seqnames = paste0("chr", genes$chromosome_name),
  #                     ranges = IRanges::IRanges(start = genes$start_position, end = genes$end_position),
  #                     gene_name = genes$hgnc_symbol, gene_biotype = genes$gene_biotype, strand = genes$strand)
  
  genes_gr <- get_gene_info_biomart()
  
  promoters_gr <- GenomicRanges::promoters(genes_gr, upstream = 1500, downstream = 500)
  
  ICR_summary <- res_pre %>% 
    dplyr::group_by(ICR_id ) %>% 
    dplyr::summarise(ICR_n = dplyr::n(), 
              ICR_chr = unique(ICR_chr), 
              ICR_start = min(ICR_start), 
              ICR_stop = max(ICR_stop), 
              mean_logFC = mean(logFC, na.rm = TRUE), 
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
    dplyr::group_by(ICR_id, ICR_n, ICR_chr, ICR_start, ICR_stop, mean_logFC, direction, n_hyper, n_hypo, perc_signif, FDR) %>% 
    dplyr::summarise(Gene = paste(na.omit(unique(Gene)), collapse = ", "), .groups = "drop") %>% 
    dplyr::arrange(FDR) %>% 
    dplyr::mutate(perc_signif = paste0(round(perc_signif, 1), "%"),
           perc_diff = ifelse(mean_logFC >= 0, paste0("+", round(mean_logFC * 100, 2), "%"),
                              paste0("-", round(abs(mean_logFC * 100), 2), "%"))) %>% 
    dplyr::relocate(perc_diff, .after = mean_logFC)
  
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
      data = ICR_annotated_plot[which(ICR_annotated_plot$FDR<sig),],aes(x=bp_cum, y = -log10(FDR), label=ICR_id),
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
    dplyr::group_by(ICR_id, ICR_n, ICR_chr, ICR_start, ICR_stop, mean_logFC, direction, n_hyper, n_hypo, perc_signif, FDR) %>%
    dplyr::summarise(Gene = paste(na.omit(unique(Gene)), collapse = ", "), .groups = "drop") %>%
    dplyr::arrange(FDR) %>%
    dplyr::mutate(
      perc_signif = paste0(round(perc_signif, 1), "%"),
      perc_diff = ifelse(mean_logFC >= 0,
                         paste0("+", round(mean_logFC * 100, 2), "%"),
                         paste0("-", round(abs(mean_logFC * 100), 2), "%"))
    ) %>%
    dplyr::relocate(perc_diff, .after = mean_logFC)
  
  
  out <- list(ICR_summary = ICR_summary, ICR_collapsed = ICR_collapsed,
              plots= list(manhplot_icr = manhplot_icr))
  # save(ICR_collapsed, file = paste0(output_dir_path, "/res_EWAS_", label, "_ICR.Rdata"))
  
  return(out)
  
  
}




#'get_gene_info_biomart
#'
#' @description Cache biomart results in path associated with this package, this avoids
#' querying biomart every time the code is run (will results with timeout/ blocks
#' from too many requests)
#' @return table of genes and associated metadata
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


