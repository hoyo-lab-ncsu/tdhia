



#' add_gene_metadata_to_icr
#'
#' @param df_res dataframe that is output from analyze_association()
#' @param sort_adj_pval sort rows by adjusted p-value
#' @export
add_gene_metadata_to_icr <- function(df_res, sort_adj_pval = TRUE) {
  
  # Filter Significance ICRs and add nearest gene
  #_______________________________________________________________________________
  
  # Significant ICRs
  df_res_sig <- filter(select(df_res, -Formula), ADJ_P_VAL <0.05)
  df_res_sig$icr_name <- stringr::str_replace(df_res_sig$icr_id, "_", " ")
  
  # Switch to pretty names for Response variables
  # df_res_sig$Response <- sapply(df_res_sig$Response, function (x)
  #     response_varnames[which(x == response_vars)])
  
  idat_dir_path <- system.file("data-raw", package = "tdhia")
  
  
  # Read list of ICRs that overlap with nearby zing finger
  df_zinc_finger <-read.csv(paste0( system.file("data-raw", package = "tdhia"), "/ICR_ZincFinger_w-in_1000bp.csv"))
  df_zinc_finger$icr_id <- as.numeric(gsub("ICR_([0-9]+).*","\\1",df_zinc_finger$icr))
  df_zinc_finger$near_zinf_finger <- rowSums(!df_zinc_finger[, 3:5] == "", na.rm = TRUE) > 0
  zinc_finger_icrs <- df_zinc_finger$icr_id[df_zinc_finger$near_zinf_finger]
  
  
  # Load annotated list of whole imprintome
  # ICR_IDs that end with "#" are high confidence
  imp_whole <- read.csv(paste0(system.file("data-raw", package = "tdhia"), "/whole_imprintome_table.csv")) 
  imp_whole$icr_id <- as.numeric(gsub("ICR_([0-9]+).*","\\1",imp_whole$ID))
  imp_whole$icr_name <- paste0("ICR ", imp_whole$icr_id)
  # Scan for previously published icrs
  high_conf_icrs <- imp_whole$icr_id[grep('#', imp_whole$ID)]
  
  # Get list of ICRs that have evidence of gametic origin for methylation
  # Includes "high confidence" (lit validated ICRs)
  imp_gamete <- read.csv(paste0(system.file("data-raw", package = "tdhia"), "/imprintome_icrs_gametic_methyl_origin.csv"))
  imp_gamete$icr_id <- as.numeric(gsub("ICR_([0-9]+).*","\\1",imp_gamete$ID))
  med_conf_icrs <- imp_gamete$icr_id #setdiff(imp_gamete$icr_id, high_conf_icrs)
  low_conf_icrs <- setdiff(imp_whole$icr_id, union(med_conf_icrs, high_conf_icrs))
  
  
  # Add zinc finger info
  df_res_sig$is_icr_zinc <- df_res_sig$icr_id %in% paste0("ICR_", zinc_finger_icrs)
  df_res_sig$is_icr_highconf <-  df_res_sig$icr_id %in% paste0("ICR_", high_conf_icrs)
  df_res_sig$is_icr_medconf <-  df_res_sig$icr_id %in% paste0("ICR_", med_conf_icrs)
  df_res_sig$is_icr_lowconf <-  df_res_sig$icr_id %in% paste0("ICR_", low_conf_icrs)
  
  
  
  # Add closest gene to each ICR
  df_res_sig <- merge(x = df_res_sig, y = imp_whole %>%
                        select(icr_name, Genomic.Coordinates, Nearest.Transcript, 
                               Distance.to.Nearest.Transcript),
                      by = "icr_name", all.x = TRUE, all.y = FALSE, sort = FALSE)
  
  if (sort_adj_pval) df_res_sig <- dplyr::arrange(df_res_sig, ADJ_P_VAL)
  
  return(df_res_sig)
}