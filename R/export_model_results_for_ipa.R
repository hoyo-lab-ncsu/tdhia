

#' export_model_results_for_ipa
#'
#'
#' IPA
#' Create New \> Upload data
#' Columns: Expr Other, Expr Fale Discovery, ENSEMBL ID
#' Only include genes below threshold in analysis file
#' Reference is IPA gene only dataset, or user dataset
#' 
#' IPA Core Analysis Dataset Mapping
#' \> column_in_csv_file: variable in IPA analysis \<\<
#'   Nearest Transcript: Gene Symbol
#'   ADJ_P_VAL: Expr False Discovery Rate (q-value)
#'   emsenbl_gene_id: Ensembl
#'   nEstimate: Expr Log Ratio
#'
#' @param df_models dataframe of imprintime associated model output.
#' @param output_dir_path path to export results from each model
#' @param imp_type specify whether the analysis done is at the icr level or 
#' cpg level.
#' @param add_ensembl_ids when trye, attempts to add corresponding ensembl IDS 
#' based ont he gene symbols.
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @export
export_model_results_for_ipa <- function(df_models, output_dir_path, imp_type = "cpg", add_ensembl_ids = F) {


df_models <- cbind(df_models, add_metadata_to_imp_sites(
  df_models$cpg_id, imp_type = "cpg",db_flag = T) %>% dplyr::select(-c("icr_id","cpg_id")))


# Significant ICRs
df_gene_unsplit <- df_models %>% 
  dplyr::select(c("model_group", "Nearest.Transcript", "Estimate" , "ADJ_P_VAL"))
df_gene_unsplit$ind = 1:nrow(df_gene_unsplit)



# Split up entries with multiple genes separated by "|"
#______________________________________________________________________________
dataframe_split_genes <- function(df) {
  temp_list = list(); k = 1
  for (n in 1:nrow(df)) {
    splits <- str_split(string = df$Nearest.Transcript[n],pattern = "\\|")
    # Add row for each entry
    for (m in seq_along(splits[[1]])) {
      temp_list[[k]] <- df[n,] %>% mutate(Nearest.Transcript = splits[[1]][m])
      k = k+1
    }
  }
  df_out <- do.call(rbind, temp_list); return(df_out) 
}
df_gene_split <- dataframe_split_genes(df_gene_unsplit)


temp_list = list(); k = 1
for (n in 1:nrow(df_gene_unsplit)) {
  splits <- str_split(string = df_gene_unsplit$Nearest.Transcript[n],pattern = "\\|")
  # Add row for each entry
  for (m in seq_along(splits[[1]])) {
    temp_list[[k]] <- df_gene_unsplit[n,] %>% mutate(Nearest.Transcript = splits[[1]][m])
    k = k+1
  }
}
df_gene_split <- do.call(rbind, temp_list) 


# Grab max magnitude model estimate for each gene (and model group)
df_genes <- df_gene_split %>% 
  group_by(model_group, Nearest.Transcript) %>% 
  summarize(Estimate = Estimate[which.max(abs(Estimate))],
            ADJ_P_VAL = min(ADJ_P_VAL), ind = mean(ind)) %>% 
  arrange(ind) %>% dplyr::select(-ind)
df_genes %>% arrange(ADJ_P_VAL)



# Attempt to find matching ensembl ids
if (add_ensembl_ids) {
  if (!requireNamespace("biomaRt", quietly = TRUE))  BiocManager::install("biomaRt", force = TRUE)
  library(biomaRt)
  
  # Connect to Ensembl (human genes)
  mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  # Query Ensembl for Ensembl ID + gene type
  mapping <- getBM(
    attributes = c("hgnc_symbol", "ensembl_gene_id", "gene_biotype"),
    filters = "hgnc_symbol",
    values = df_genes$Nearest.Transcript,
    mart = mart
  )
  # Use first entry only from biomart
  mapping_unique <- mapping %>%  dplyr::group_by(hgnc_symbol) %>%
    dplyr::slice(1) %>%  dplyr::ungroup()
  
  
  # Add biomart matches to gene list dataframe
  df_genes <- df_genes %>% left_join(
    mapping_unique, by = join_by(Nearest.Transcript == hgnc_symbol ), keep = F)
  
  
  # Attempt to fill in genes with annotation DBI
  remaining_genes_bv = is.na(df_genes$ensembl_gene_id)
  ensemble_ids = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                       keys=df_genes$Nearest.Transcript[remaining_genes_bv],
                                       column="ENSEMBL",
                                       keytype="SYMBOL",
                                       multiVals="first")
  df_genes$ensembl_gene_id[remaining_genes_bv] <- ensemble_ids
  df_genes$negEstimate = -df_genes$Estimate
  # df_genes$negInvEstimate = -(1/df_genes$Estimate)
}


exported_datasets <- list()
model_names <- unique(df_genes$model_group)
for (n in seq_along(model_names)) {
  model_path = paste0(output_dir_path, "/", model_names[n],"-group__gene_pathway_analysis_input_all_icrs.csv")
  cat(sprintf("> Exporting: %s\n", model_path))
  exported_datasets[[n]] <- df_genes %>% filter(model_group==model_names[n]) %>% dplyr::arrange(ADJ_P_VAL)
  write.csv(x = exported_datasets[[n]], file = model_path)
  
}



# Get list of all genes on imprintome
# df_gene_unsplit_imprintome <- tdhia::imprintome_icr_nearest_transcripts %>% dplyr::select(ID, Nearest.Transcript)
# df_gene_split_imprintome <- dataframe_split_genes(df_gene_unsplit_imprintome)
# 
# write.csv(x = df_gene_split_imprintome, file =  paste0(output_dir_path, "/gene_pathway_analysis_all_genes_reference.csv"))

return(list(exported_datasets = exported_datasets))


}


