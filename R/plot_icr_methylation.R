


#' plot_icr_dotplot
#' @description produces a simple dot pot of beta values of CpG sites within a 
#' specified ICR, seperating patients into 2 groups.
#' 
#' @param mat_cpg_beta a matrix of beta values, cpg sites (rows) x patients (col)
#' @param sig_cpgs vector of cpg site IDs that are significant (can be across 
#' all ICRs and not just ICR being plotted).
#' @param df_patient_groups a dataframe with the sample_colname 
#' (default: patient_id) and a "group" column that specifies whether each patient 
#' is in the low (1) or high (2) group.
#' @param icr_id the ID of the specific ICR to be plotted of the form ICR_#, i.e., ICR_10.
#' @param xlab_txt xlabel for plot 
#' @param plot_height_width vector of two numbers in inches, for plot width and 
#' height, as defined in cowplot::save_plot().
#' @param max_sig_hwindow dictates the window of CPG sites to be included,
#'   1) NULL: plot all cpg sites (default).
#'   2) Vector of 2 numerics: the lower and upper index of CpG sites to be 
#'   included. If positive, the index is from the start. If negative, index is 
#'   from the end. 
#'   3) Single numeric: how much padding is added around significant CpG sites.0
#' @param output_path full file system path so save plot to.
#' @param db_flag when true, saves environment variables to disk to a file with 
#' same name as function.
#' @param filter_na_group boolean, when TRUE (default), removes NA group from plots.
#' @param legend.position specify position of legend as specified in ggplot theme
#'  (default = "none").
#' @param manual_cpg_index todo
#' @param sample_colname name of coolumn that refers to sample_id found with 
#' df_patient_groups, and also matc the column names in mat_cpg_beta.
#' 
#' @returns a ggplot function handle to the plot.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export
plot_icr_dotplot <- function(mat_cpg_beta, sig_cpgs = NA, df_patient_groups, icr_id, xlab_txt = "", 
                                 plot_height_width = c(5,3), output_path, max_sig_hwindow = NULL, db_flag = F,
                                 filter_na_group = T, legend.position = "none", ytext = "Mean Beta Value", 
                                 sample_colname = "patient_id", overwrite_plot = F) {
  # Create output folder
  dir.create(path = output_path, recursive = TRUE, showWarnings = FALSE)
  if(db_flag) save(list = ls(all.names = TRUE), file = "plot_icr_methylation.RData")
  # load(file = "plot_icr_methylation.RData")
  
  # Get list of CpGs for specified icrs
  df = left_join(x = data.frame(cpg_id = rownames(mat_cpg_beta)),
            y = tdhia::manifest_v1A2_design_scores %>% dplyr::select("cpg_id", "icr_id","MAPINFO") %>% distinct(),
            by = join_by("cpg_id"), keep = FALSE, na_matches = "never", 
            relationship = "one-to-one")
  
  # Subset the cpg_beta matrix
  sub_mat_cpg_beta = mat_cpg_beta[df$icr_id == icr_id,]
  
  # Reorder rows to genomic location
  sub_mat_cpg_beta <- sub_mat_cpg_beta %>% arrange(df[df$icr_id == icr_id,]$MAPINFO)
  # ordered_cpg_ids <- df[df$icr_id == icr_id,] %>% arrange(MAPINFO) %>% pull(cpg_id)
  
  # Subset the cpg sites to those around significant cpgs/ specified by user
  if (length(max_sig_hwindow)==2) {
    # If indices are negative, assume index is from end, convert to positive from start
    if (max_sig_hwindow[1] < 0) max_sig_hwindow[1] = nrow(sub_mat_cpg_beta) + max_sig_hwindow[1]
    if (max_sig_hwindow[2] < 0) max_sig_hwindow[2] = nrow(sub_mat_cpg_beta) + max_sig_hwindow[2]
    
    
    sub_mat_cpg_beta <- sub_mat_cpg_beta[max_sig_hwindow[1]:max_sig_hwindow[2],]
  } else if (!is.null(max_sig_hwindow) && (length(max_sig_hwindow)>0) &&
             (nrow(sub_mat_cpg_beta) > 2*max_sig_hwindow+1)) {
    
    sig_inds <- which(rownames(sub_mat_cpg_beta) %in% sig_cpgs)
    
    min_ind <- max(c(min(sig_inds) - max_sig_hwindow, 1))
    max_ind <- min(c(max(sig_inds) + max_sig_hwindow, nrow(sub_mat_cpg_beta)))
    
    sub_mat_cpg_beta <- sub_mat_cpg_beta[min_ind:max_ind,]
    
  }

  # Ordered list of cpg sites (for factor)
  # rownames(sub_mat_cpg_beta)
  # Convert data to long format
  df_long_cpg_beta <- sub_mat_cpg_beta %>% rownames_to_column("cpg_id") %>% 
    pivot_longer(cols = -"cpg_id", names_to = sample_colname)
  
  df_long_cpg_beta <- left_join(x = df_long_cpg_beta, y= df_patient_groups, by = join_by({{sample_colname}}),
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
  if (filter_na_group) df_summary <- df_summary %>% filter(!is.na(group))
  
  # Get ICR metadata (closest genes, zinc finger)
  icr_metadata <- add_metadata_to_imp_sites(icr_id, imp_type = "icr")
  zinc_finger_str = c("-","+")[(as.numeric(icr_metadata$is_icr_zinc)+1)]
  icf_conf_str = c("High","Medium", "Low")[icr_metadata$icr_conf]
  
  
  exact_xlim = range(df_summary$beta_mean)
  padded_xlim  = c(exact_xlim[1] - 0.1*diff(exact_xlim),exact_xlim[2] + 0.1*diff(exact_xlim))
  
 
  # Plot methylation across cpg sites in ICR
  gg <- ggplot(data = df_summary, aes(x = cpg_id, y = beta_mean)) +
    geom_point(aes(color = group), size = 1.5, alpha = 0.5) + 
    geom_text(aes(label = ifelse(df_summary$cpg_id %in% sig_cpgs, "*", "")), 
              y = padded_xlim[2], size = 4, vjust=1) + 
    # geom_boxplot(aes(color = group), position = position_dodge(width=0.7), width = 0.7, linewidth = .2) +
    scale_color_manual(values=c("blue", "red")) + #, labels= c("Low","High"), guide = "none") + 
    scale_fill_manual(values=c("grey90", "white"), guide = "none") +
    # scale_x_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = (c(min(df_summary$cpg_id_rank), max(df_summary$cpg_id_rank))), ylim = padded_xlim) + #ylim=c(0,1)
      xlab(xlab_txt) + ylab(ytext) + 
    # geom_ribbon(aes(fill = ))+
    # ggtitle(sprintf("%s (%sZF, %s): %s", icr_id, zinc_finger_str,icf_conf_str, 
    #           icr_metadata$Nearest.Transcript)) +
    ggtitle(sprintf("%s: %s", icr_id, icr_metadata$Nearest.Transcript)) +
    theme_classic(base_size = 7) + theme(axis.text.x = element_text(
      angle = 45,vjust = 1, hjust = 1), plot.title = element_text(size = 7),
      legend.position = legend.position) 
  
  
  plot_path <- paste0(output_path, "/", 
                      sprintf("Beta_%sZF_%s_%s", zinc_finger_str, icf_conf_str, icr_id), ".jpg")
  if (overwrite_plot) {
    # Print summaries of data to command line as well
    print(gg)
    print(table(df_patient_groups$group))
    save_plot(filename = plot_path, plot = gg,base_height = 2, base_width = 2.5)
  }

  
  # Export
  return(list(plot = gg, df_summary = df_summary, cpg_beta_plotted = sub_mat_cpg_beta))
}







#' plot_icr_diffbar
#' @description produces a simple dot pot of beta values of CpG sites within a 
#' specified ICR, seperating patients into 2 groups.
#' 
#' Note: cpg sites may appear out of order, but that is because they are sorted 
#' by their genomic position (MAPINFO) and not label
#' 
#' @param mat_cpg_beta a matrix of beta values, cpg sites (rows) x patients (col)
#' @param sig_cpgs vector of cpg site IDs that are significant (can be across 
#' all ICRs and not just ICR being plotted).
#' @param df_patient_groups a dataframe with the  
#' 1) a sample colum  name (default: patient_id, specified by sample_colname) 
#' 2) diff_group: column that specifies whether each patient 
#' is in the low (1) or high (2) group.
#' 3) subset_group: column that specifies different subsets of patient population
#'  (combined will be plotted also).
#' @param icr_id the ID of the specific ICR to be plotted of the form ICR_#, i.e., ICR_10.
#' @param xlab_txt xlabel for plot 
#' @param plot_height_width vector of two numbers in inches, for plot width and 
#' height, as defined in cowplot::save_plot().
#' @param max_sig_hwindow dictates the window of CPG sites to be included,
#'   1) NULL: plot all cpg sites (default).
#'   2) Vector of 2 numerics: the lower and upper index of CpG sites to be 
#'   included. If positive, the index is from the start. If negative, index is 
#'   from the end. 
#'   3) Single numeric: how much padding is added around significant CpG sites.0
#' @param output_path full file system path so save plot to.
#' @param db_flag when true, saves environment variables to disk to a file with 
#' same name as function.
#' @param filter_na_group boolean, when TRUE (default), removes NA group from plots.
#' @param legend.position specify position of legend as specified in ggplot theme
#'  (default = "none").
#' @param manual_cpg_index todo
#' @param sample_colname name of coolumn that refers to sample_id found with 
#' df_patient_groups, and also matc the column names in mat_cpg_beta.
#' 
#' @returns a ggplot function handle to the plot.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export
plot_icr_diffbar <- function(mat_cpg_beta, sig_cpgs = NA, df_patient_groups, icr_id, xlab_txt = "", 
                             plot_height_width = c(3,2), output_path, max_sig_hwindow = NULL, db_flag = F,
                             filter_na_group = T, legend.position = "none", 
                             sample_colname = "patient_id", overwrite_plot = T) {
  # Create output folder
  dir.create(path = output_path, recursive = TRUE, showWarnings = FALSE)
  if(db_flag) save(list = ls(all.names = TRUE), file = "plot_icr_diffbar.RData")
  # load(file = "plot_icr_diffbar.RData")
  
  # Get list of CpGs for specified icrs
  df = left_join(x = data.frame(cpg_id = rownames(mat_cpg_beta)),
                 y = tdhia::manifest_v1A2_design_scores %>% dplyr::select("cpg_id", "icr_id","MAPINFO") %>% distinct(),
                 by = join_by("cpg_id"), keep = FALSE, na_matches = "never", 
                 relationship = "one-to-one")
  
  # Subset the cpg_beta matrix
  sub_mat_cpg_beta = mat_cpg_beta[df$icr_id == icr_id,]
  
  # Reorder rows to genomic location, assuming same chromosomes
  sub_mat_cpg_beta <- sub_mat_cpg_beta %>% arrange(df[df$icr_id == icr_id,]$MAPINFO)
  # ordered_cpg_ids <- df[df$icr_id == icr_id,] %>% arrange(MAPINFO) %>% pull(cpg_id)
  
  # Subset the cpg sites to those around significant cpgs/ specified by user
  if (length(max_sig_hwindow)==2) {
    # If indices are negative, assume index is from end, convert to positive from start
    if (max_sig_hwindow[1] < 0) max_sig_hwindow[1] = nrow(sub_mat_cpg_beta) + max_sig_hwindow[1]
    if (max_sig_hwindow[2] < 0) max_sig_hwindow[2] = nrow(sub_mat_cpg_beta) + max_sig_hwindow[2]
    
    
    sub_mat_cpg_beta <- sub_mat_cpg_beta[max_sig_hwindow[1]:max_sig_hwindow[2],]
  } else if (!is.null(max_sig_hwindow) && (length(max_sig_hwindow)>0) &&
             (nrow(sub_mat_cpg_beta) > 2*max_sig_hwindow+1)) {
    
    sig_inds <- which(rownames(sub_mat_cpg_beta) %in% sig_cpgs)
    
    min_ind <- max(c(min(sig_inds) - max_sig_hwindow, 1))
    max_ind <- min(c(max(sig_inds) + max_sig_hwindow, nrow(sub_mat_cpg_beta)))
    
    sub_mat_cpg_beta <- sub_mat_cpg_beta[min_ind:max_ind,]
    
  }

  # Ordered list of cpg sites (for factor)
  # rownames(sub_mat_cpg_beta)
  # Convert data to long format
  df_long_cpg_beta <- sub_mat_cpg_beta %>% rownames_to_column("cpg_id") %>% 
    pivot_longer(cols = -"cpg_id", names_to = sample_colname)
  
  df_long_cpg_beta <- left_join(x = df_long_cpg_beta, y= df_patient_groups, by = join_by({{sample_colname}}),
                                keep = FALSE, na_matches = "never", relationship = "many-to-one")
  
  df_long_cpg_beta$cpg_id <- factor(df_long_cpg_beta$cpg_id, levels = rev(rownames(sub_mat_cpg_beta)), ordered=TRUE)
  df_long_cpg_beta$diff_group <- factor(df_long_cpg_beta$diff_group, levels = unique(df_long_cpg_beta$diff_group) %>% sort(), ordered=TRUE)
  df_long_cpg_beta$subset_group <- factor(df_long_cpg_beta$subset_group, levels = unique(df_long_cpg_beta$subset_group) %>% sort(), ordered=TRUE)
  
  
  # Calculate beta difference across all patients
  df_summary_all <- df_long_cpg_beta %>% group_by(cpg_id) %>% 
    summarize(
      subset_group = "All",
      beta_mean_diff = mean(value[diff_group==2], na.rm = T) - mean(value[diff_group==1], na.rm = T),
      n1 = sum(diff_group==1), n2 = sum(diff_group==2),
      beta_sd_diff = sqrt(sd(value[diff_group==1], na.rm = T)^2 +
                            sd(value[diff_group==2], na.rm = T)^2),
      beta_sem_diff = sqrt(sd(value[diff_group==1], na.rm = T)^2/n1 +
                             sd(value[diff_group==2], na.rm = T)^2/n2))
  df_summary_all$cpg_id_rank = as.numeric(df_summary_all$cpg_id)
  df_summary_all$xmin = df_summary_all$cpg_id_rank -0.5
  df_summary_all$xmax = df_summary_all$cpg_id_rank +0.5
  df_summary_all$back_fill = df_summary_all$cpg_id_rank %% 2 == 0
  df_summary_all$cpg_sig = df_summary_all$cpg_id %in% sig_cpgs
  
  # Get difference in beta for each of the sample subgroups
  df_summary_sub <- df_long_cpg_beta %>% group_by(cpg_id, subset_group) %>% 
    summarize(
      beta_mean_diff = mean(value[diff_group==2], na.rm = T) - mean(value[diff_group==1], na.rm = T),
      n1 = sum(diff_group==1), n2 = sum(diff_group==2),
      beta_sd_diff = sqrt(sd(value[diff_group==1], na.rm = T)^2 +
                            sd(value[diff_group==2], na.rm = T)^2),
      beta_sem_diff = sqrt(sd(value[diff_group==1], na.rm = T)^2/n1 +
                             sd(value[diff_group==2], na.rm = T)^2/n2))
  
  temp <- df_summary_all %>% dplyr::select(c("cpg_id_rank", "xmin", "xmax", "back_fill", "cpg_sig"))
  df_summary_sub <- cbind(df_summary_sub, rbind(temp,temp))
  # Bind summary stat of all group and subset groups
  df_summary <- rbind(df_summary_all, df_summary_sub)
  # Enforce factor level order (reverse because y-axis inverted in plotting)
  df_summary$subset_group <- factor(
    df_summary$subset_group, levels = rev(c("All", levels(df_long_cpg_beta$subset_group))), ordered = TRUE)

  
  # Get ICR metadata (closest genes, zinc finger)
  icr_metadata <- add_metadata_to_imp_sites(icr_id, imp_type = "icr")
  zinc_finger_str = c("-","+")[(as.numeric(icr_metadata$is_icr_zinc)+1)]
  icf_conf_str = c("High","Medium", "Low")[icr_metadata$icr_conf]

  exact_xlim = range(df_summary$beta_mean_diff)
  padded_xlim  = c(exact_xlim[1] - 0.5*diff(exact_xlim),exact_xlim[2] + 0.5*diff(exact_xlim))
  

  # Plot methylation across cpg sites in ICR
  gg <- ggplot(data = df_summary, aes(y = cpg_id, x = beta_mean_diff, )) +
    geom_rect(data = df_summary_all,  aes(ymin = xmin, ymax = xmax, xmin = -Inf, xmax = +Inf), 
               fill = ifelse(df_summary_all$back_fill, "grey92", "white"),
               color = ifelse(df_summary_all$cpg_sig, "black", NA), linewidth = 0.25) +
    geom_col(aes(fill = subset_group),alpha = 1, position = "dodge", width = 1) + 
    scale_fill_manual(values = c("All" = "black", "1" = "#f03b20", "2" = "#67a9cf"))+
    coord_cartesian(ylim = c(min(df_summary$cpg_id_rank)-0.5,
                             max(df_summary$cpg_id_rank)+0.5), xlim = padded_xlim, expand = c(0,0)) +
    geom_vline(xintercept = 0, color = "black")+
    ylab(xlab_txt) + xlab("Mean Beta Value") + 
    ggtitle(sprintf("%s: %s", icr_id, icr_metadata$Nearest.Transcript)) +
    theme_classic(base_size = 7) + 
    theme(axis.text.x = element_text(vjust = 0.5, hjust = 1),
          axis.text.y = element_text(colour = ifelse( df_summary$cpg_sig, "black", "grey50"),
                                     face = "bold"),
      plot.title = element_text(size = 7),
      legend.position = legend.position)
  gg
  
  
  plot_path <- paste0(output_path, "/", 
                      sprintf("Beta_%sZF_%s_%s", zinc_finger_str, icf_conf_str, icr_id), ".jpg")
  if (overwrite_plot | !file.exists(plot_path)) {
    # Print summaries of data to command line as well
    print(gg)
    print(table(df_summary$diff_group))
    save_plot(filename = plot_path, plot = gg, base_height = plot_height_width[1], base_width = plot_height_width[1])
  }

  
  # Export
  return(list(plot = gg, df_summary = df_summary, cpg_beta_plotted = sub_mat_cpg_beta))
}
