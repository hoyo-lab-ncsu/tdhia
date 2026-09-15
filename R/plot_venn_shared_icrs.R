
#' Plot Shared ICRs Across Groups
#'
#' Builds a Venn diagram from unique ICR identifiers in each group, with
#' groups ordered by decreasing number of distinct ICRs.
#'
#' @param df_imp_sig Data frame containing icr_id and the grouping column.
#'   Supply already-selected significant sites; no significance filter is
#'   applied by this function.
#' @param colname Character string naming the grouping column.
#' @param group_names Optional vector of group values to include. NULL uses
#'   unique values of the grouping column. These values select groups;
#'   they are not replacement display labels.
#' @param output_dir_path Existing output directory for the optional PNG.
#'   Defaults to the working directory.
#' @param name_suffix Optional filename suffix. NULL disables saving;
#'   any non-NULL value, including an empty string, enables PNG output.
#' @param db_flag Logical; save the initial environment to
#'   plot_venn_shared_icrs.RData in the working directory.
#'
#' @details
#' When saving is enabled, the filename is icr_venn_diagram followed by
#' name_suffix and .png. The plot is saved at a base width of four inches
#' and height of two inches. The directory is not created by this function.
#'
#' @return A named list with enn_input, the named list of unique ICR vectors
#'   used for the diagram, and figure, the plot object returned by
#'   ggVennDiagram::ggVennDiagram() with additional theme settings.
#'   The field name enn_input reflects the current implementation.
#' @export
plot_venn_shared_icrs <- function(df_imp_sig, colname = "Response", group_names = NULL, 
                             output_dir_path = getwd(), name_suffix = NULL, db_flag = FALSE) {
  
  if (db_flag) {save(list = ls(all.names = TRUE), file = "plot_venn_shared_icrs.RData")}
  # load(file = "plot_venn_shared_icrs.RData")
  
  # If not specified groups for the venn diagram are found in the group_names column
  if (is.null(group_names)) group_names = unique(df_imp_sig[[colname]])
  
  # For each group included in venn diagram, get list of all ICRs in that group
  df_members = list()
  for (n in seq_along(group_names)) {
    df_members[[group_names[n]]] <- unique(dplyr::filter(df_imp_sig, .data[[colname]] == group_names[n])$icr_id)
  }
  # Sort groups by total number of members
  df_members = df_members[ order(sapply(df_members, function(x) length(x)),decreasing = TRUE)]
 
  # Input to venn diagram package
  # browser()
  # png(filename=paste0(output_dir_path, "/icr_venn_diagram",name_suffix,".png"))
  ggv = ggVennDiagram::ggVennDiagram(df_members, label = "count", category.names = 
                        names(df_members),label_alpha = 1,set_color = "black", label_color = "black",
                      set_size = 4, label_size = 3.5, order.intersect.by = "size") + 
    scale_x_continuous(expand = expansion(mult = .2)) +
    theme(legend.text=element_text(size=8), legend.title = element_text(size=8),
          legend.position = "none")
  if (!is.null(name_suffix)) {
    cowplot::save_plot(filename = paste0(output_dir_path, "/icr_venn_diagram",name_suffix,".png"),
              plot = ggv,base_height = 2,base_width = 4)
  } 
  
  # if (!is.null(name_suffix)) dev.off()
  ggv
  
  return(list(enn_input = df_members, figure = ggv)) #v, 
}
