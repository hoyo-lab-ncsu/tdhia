
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
                             output_dir_path = getwd(), name_suffix = NULL) {
  
  # browser()
  # If not specified groups for the venn diagram are found in the group_names column
  if (is.null(group_names)) group_names = unique(df_imp_sig[[colname]])
  
  # For each group included in venn diagram, get list of all ICRs in that group
  df_members = list()
  for (n in seq_along(group_names)) {
    df_members[[group_names[n]]] <- unique(filter(df_imp_sig, .data[[colname]] == group_names[n])$icr_id)
  }
  
  # Input to venn diagram package
  
  if (!is.null(name_suffix)) png(filename=paste0(output_dir_path, "/icr_venn_diagram",name_suffix,".png"))
  ggv = ggVennDiagram(df_members, label = "count", category.names = 
                  names(df_members), set_size = 8, label_size = 8,) + 
    theme(legend.text=element_text(size=15), legend.title = element_text(size=15))
  if (!is.null(name_suffix)) dev.off()
  ggv
  
  return(list(venn_input = df_members, figure = ggv))
}