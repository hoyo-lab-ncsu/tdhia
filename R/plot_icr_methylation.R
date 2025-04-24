


#' plot_icr_methylation
#'
#'  position.
#' @param cpg_beta_mat a matrix of beta values, cpg sites (rows) x patients (col)
#' @param patient_groups n
#' @param icr_id v
#' @plot_height_width v
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @returns a ggplot function handle to the plot.
#'
#' @export
plot_icr_methylation <- function(mat_cpg_beta, df_patient_groups, icr_id, 
                                 plot_height_width = c(5,3)) {
  
  # Get list of CpGs for specified icr
  df = left_join(x = data.frame(cpg_id = rownames(mat_cpg_beta)),
            y = tdhia::manifest_v1A2_design_scores %>% select("cpg_id", "icr_id","MAPINFO") %>% distinct(),
            by = join_by("cpg_id"), keep = FALSE, na_matches = "never", 
            relationship = "one-to-one")
  

  # Subset the cpg_beta matrix
  sub_mat_cpg_beta = mat_cpg_beta[df$icr_id == icr_id,]
  
  # Reorder rows to genomic location
  sub_mat_cpg_beta <- sub_mat_cpg_beta %>% arrange(df[df$icr_id == icr_id,]$MAPINFO)
  
  ordered_cpg_ids <- df[df$icr_id == icr_id,] %>% arrange(MAPINFO) %>% pull(cpg_id)
  
  
  
  
  
  # Ordered list of cpg sites (for factor)
  # rownames(sub_mat_cpg_beta)
  # Convert data to long format
  df_long_cpg_beta <- sub_mat_cpg_beta %>% rownames_to_column("cpg_id") %>% 
    pivot_longer(cols = -"cpg_id", names_to = "patient_id")
  
  df_long_cpg_beta <- left_join(x = df_long_cpg_beta, y= df_patient_groups, by = join_by("patient_id"),
            keep = FALSE, na_matches = "never", relationship = "many-to-one")
  
  df_long_cpg_beta$cpg_id <- factor(df_long_cpg_beta$cpg_id, levels = ordered_cpg_ids, ordered=TRUE)
  df_long_cpg_beta$group <- factor(df_long_cpg_beta$group, levels = unique(df_long_cpg_beta$group) %>% sort(), ordered=TRUE)
  
  df_summary <- df_long_cpg_beta %>% group_by(cpg_id, group) %>% 
    summarize(beta_mean = mean(value, na.rm = T), beta_sd = sd(value, na.rm = T),
              beta_sd = sd(value, na.rm = T), beta_sem = sd(value, na.rm = T)/sqrt(n()))
  
  # Plot methylation across cpg sites in ICR
  ggplot(data = df_summary, aes(x = cpg_id, y = beta_mean, color = group)) +
    geom_errorbar(aes(ymin = beta_mean - beta_sem, ymax = beta_mean + beta_sem),
                  position = position_dodge(width=0.6), width = 0.4) +
    geom_boxplot(position = position_dodge(width=0.6), width = 0.5) + 
    scale_color_manual(values=c("red", "blue")) + 
    coord_cartesian(ylim=c(0,1)) +xlab("CpG site") + ylab("Beta Value") + 
    guides(color=guide_legend(title="Med Diet Score")) +
    theme_classic()
  
  # Export
}