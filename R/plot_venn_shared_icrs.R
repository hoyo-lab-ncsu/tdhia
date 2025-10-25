
#' venn_shared_icrs
#' @description produces venn diagram plot of shared ICRs
#'
#' @param df_summary a dataframe of significant cpg sites for each icr. 
#'  Required to have the (colname) specified and icr_id column.
#' @param colname column name to speficy venn diagram groups
#' @param group_names manual override of groups for venn diagram.
#' @param output_dir_path path to output folder if export is specified.
#' @param name_suffix estra misc. suffic to saved filename if export is specified.
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
    df_members[[group_names[n]]] <- unique(filter(df_imp_sig, .data[[colname]] == group_names[n])$icr_id)
  }
  # Sort groups by total number of members
  df_members = df_members[ order(sapply(df_members, function(x) length(x)),decreasing = TRUE)]
 
  # Input to venn diagram package
  # browser()
  # png(filename=paste0(output_dir_path, "/icr_venn_diagram",name_suffix,".png"))
  ggv = ggVennDiagram(df_members, label = "count", category.names = 
                        names(df_members),label_alpha = 1,set_color = "black", label_color = "white",
                      set_size = 4, label_size = 3.5, order.intersect.by = "size") + 
    scale_x_continuous(expand = expansion(mult = .2)) +
    theme(legend.text=element_text(size=8), legend.title = element_text(size=8),
          legend.position = "none")
  if (!is.null(name_suffix)) {
    save_plot(filename = paste0(output_dir_path, "/icr_venn_diagram",name_suffix,".png"),
              plot = ggv,base_height = 2,base_width = 4)
  } 
  
  # if (!is.null(name_suffix)) dev.off()
  ggv
  
  return(list(enn_input = df_members, figure = ggv)) #v, 
}