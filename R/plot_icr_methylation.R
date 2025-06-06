


#' plot_icr_methylation
#'
#'  position.
#' @param cpg_beta_mat a matrix of beta values, cpg sites (rows) x patients (col)
#' @param patient_groups n
#' @param icr_id v
#' @param plot_height_width v
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @returns a ggplot function handle to the plot.
#'
#' @export
plot_icr_methylation <- function(mat_cpg_beta, sig_cpgs = NA, df_patient_groups, icr_id, xlab_txt = "", 
                                 plot_height_width = c(5,3), output_path, max_sig_hwindow = 9) {
  
  
  # Create output folder
  dir.create(path = output_path, recursive = TRUE, showWarnings = FALSE)
  
  # Get list of CpGs for specified icr
  df = left_join(x = data.frame(cpg_id = rownames(mat_cpg_beta)),
            y = tdhia::manifest_v1A2_design_scores %>% dplyr::select("cpg_id", "icr_id","MAPINFO") %>% distinct(),
            by = join_by("cpg_id"), keep = FALSE, na_matches = "never", 
            relationship = "one-to-one")
  
  # Subset the cpg_beta matrix
  sub_mat_cpg_beta = mat_cpg_beta[df$icr_id == icr_id,]
  
  # Reorder rows to genomic location
  sub_mat_cpg_beta <- sub_mat_cpg_beta %>% arrange(df[df$icr_id == icr_id,]$MAPINFO)
  # ordered_cpg_ids <- df[df$icr_id == icr_id,] %>% arrange(MAPINFO) %>% pull(cpg_id)
  
  # Subset the cpg sites to those around significant cpgs
  if (length(max_sig_hwindow)==2) {
    if (max_sig_hwindow[1] < 0) max_sig_hwindow[1] = nrow(sub_mat_cpg_beta) + max_sig_hwindow[1]
    if (max_sig_hwindow[2] < 0) max_sig_hwindow[2] = nrow(sub_mat_cpg_beta) + max_sig_hwindow[2]
    
    sub_mat_cpg_beta <- sub_mat_cpg_beta[max_sig_hwindow[1]:max_sig_hwindow[2],]
  } else if (!is.na(max_sig_hwindow) & nrow(sub_mat_cpg_beta) > 2*max_sig_hwindow+1) {
    
    sig_inds <- which(rownames(sub_mat_cpg_beta) %in% sig_cpgs)
    
    min_ind <- max(c(min(sig_inds) - max_sig_hwindow, 1))
    max_ind <- min(c(max(sig_inds) + max_sig_hwindow, nrow(sub_mat_cpg_beta)))
    
    sub_mat_cpg_beta <- sub_mat_cpg_beta[min_ind:max_ind,]
    
  }

  
  # Ordered list of cpg sites (for factor)
  # rownames(sub_mat_cpg_beta)
  # Convert data to long format
  df_long_cpg_beta <- sub_mat_cpg_beta %>% rownames_to_column("cpg_id") %>% 
    pivot_longer(cols = -"cpg_id", names_to = "patient_id")
  
  df_long_cpg_beta <- left_join(x = df_long_cpg_beta, y= df_patient_groups, by = join_by("patient_id"),
            keep = FALSE, na_matches = "never", relationship = "many-to-one")
  
  df_long_cpg_beta$cpg_id <- factor(df_long_cpg_beta$cpg_id, levels = rownames(sub_mat_cpg_beta), ordered=TRUE)
  df_long_cpg_beta$group <- factor(df_long_cpg_beta$group, levels = unique(df_long_cpg_beta$group) %>% sort(), ordered=TRUE)
  
  df_summary <- df_long_cpg_beta %>% group_by(cpg_id, group) %>% 
    summarize(beta_mean = mean(value, na.rm = T), beta_sd = sd(value, na.rm = T),
              beta_sd = sd(value, na.rm = T), beta_sem = sd(value, na.rm = T)/sqrt(n()))
  df_summary$cpg_id_rank = as.numeric(df_summary$cpg_id)
  df_summary$xmin = df_summary$cpg_id_rank -0.5
  df_summary$xmax = df_summary$cpg_id_rank +0.5
  df_summary$back_fill = df_summary$cpg_id_rank %% 2 == 0
  
  
  # Get ICR metadata (closest genes, zinc finger)
  icr_metadata <- add_metadata_to_imp_sites(icr_id, imp_type = "icr")
  zinc_finger_str = c("-","+")[(as.numeric(icr_metadata$is_icr_zinc)+1)]
  icf_conf_str = c("High","Medium", "Low")[icr_metadata$icr_conf]
  
  
  exact_ylim = range(df_summary$beta_mean)
  padded_ylim  = c(exact_ylim[1] - 0.1*diff(exact_ylim),exact_ylim[2] + 0.1*diff(exact_ylim))
  
  # Plot methylation across cpg sites in ICR
  gg <- ggplot(data = df_summary, aes(x = cpg_id, y = beta_mean)) +
    # geom_rect(aes(fill = back_fill, xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf)) +
    # geom_linerange (aes(x = cpg_id, ymin = beta_mean - beta_sem,
    #                  ymax = beta_mean + beta_sem, color = group),
    #              linewidth = 5, alpha = 0.5, position = position_dodge(width=0.7))+
    geom_point(aes(color = group), size = 1.5, alpha = 0.5) + 
    # geom_errorbar(aes(ymin = beta_mean - beta_sem, ymax = beta_mean + beta_sem, color = group),
                  # position = position_dodge(width=0.7), width = 0.6, linewidth = .2) +
    geom_text(aes(label = ifelse(df_summary$cpg_id %in% sig_cpgs, "*", "")), 
              y = padded_ylim[2], size = 4, vjust=1) + 
    # geom_boxplot(aes(color = group), position = position_dodge(width=0.7), width = 0.7, linewidth = .2) +
    scale_color_manual(values=c("red", "blue")) + #, labels= c("Low","High"), guide = "none") + 
    scale_fill_manual(values=c("grey90", "white"), guide = "none") +
    # scale_x_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = (c(min(df_summary$cpg_id_rank), max(df_summary$cpg_id_rank))), ylim = padded_ylim) + #ylim=c(0,1)
      xlab(xlab_txt) + ylab("Mean Beta Value") + 
    # geom_ribbon(aes(fill = ))+
    # ggtitle(sprintf("%s (%sZF, %s): %s", icr_id, zinc_finger_str,icf_conf_str, 
    #           icr_metadata$Nearest.Transcript)) +
    ggtitle(sprintf("%s: %s", icr_id, icr_metadata$Nearest.Transcript)) +
    theme_classic(base_size = 7) + theme(axis.text.x = element_text(
      angle = 45,vjust = 1, hjust = 1), plot.title = element_text(size = 7)) 
  print(gg)
  save_plot(filename = paste0(output_path, "/", 
                              sprintf("%s_%s_%sZF", icr_id, icf_conf_str, zinc_finger_str), ".jpg"),
            plot = gg,base_height = 2, base_width = 2.5)
  
  
  
  ifelse(levels(df_summary$cpg_id) %in% NA, "purple", "black")
  # Export
}