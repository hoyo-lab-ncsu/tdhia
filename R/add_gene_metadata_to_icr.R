



#' add_metadata_to_imp_sites
#'
#' @param imp_ids list of ids for getting metadata
#' @param imp_type specify whether ids are cpg ids or icr ids
#' @param db_flag boolean, when true, save environemnt variables to disk.
#' @importFrom rlang .data
#' @export
add_metadata_to_imp_sites <- function(imp_ids, imp_type = "cpg", db_flag = T) {

  if (db_flag) {save(list = ls(all.names = TRUE), file = "add_metadata_to_imp_sites.RData")}
  # load(file = "add_metadata_to_imp_sites.RData")
  
  if (imp_type == "cpg") {
     df <- data.frame(cpg_id = imp_ids)    
    df$icr_id <- sapply(df$cpg_id, function(x)
      tdhia::mapping_cpg_icr_ids$ICR_id[x==tdhia::mapping_cpg_icr_ids$CpG_id][1])
  } else {
    df <- data.frame(icr_id = imp_ids)
  }
  # df$icr_id <- as.numeric(str_replace(df$icr_name, "^ICR ",""))
  
  # Filter Significance ICRs and add nearest gene
  #_______________________________________________________________________________


  # Read list of ICRs that overlap with nearby zing finger
  df_zf <- tdhia::imprintome_icr_zinc_finger
  
  df_zinc_finger <- data.frame(icr_id = paste0("ICR_", gsub("ICR_([0-9]+).*","\\1",df_zf$icr)))
  df_zinc_finger$ZFP57 <- rowSums(df_zf[,3-7] == "ZFP57") > 0
  df_zinc_finger$ZFP445 <- rowSums(df_zf[,3-7] == "ZFP445") > 0
  df_zinc_finger$near_zinf_finger <-   df_zinc_finger$ZFP57 |   df_zinc_finger$ZFP445
  # zinc_finger_icrs <- df_zinc_finger$icr_id[df_zinc_finger$near_zinf_finger]


  # Load annotated list of whole imprintome
  # ICR_IDs that end with "#" are high confidence
  imp_whole <- tdhia::imprintome_icr_nearest_transcripts
  imp_whole$icr_id <- paste0("ICR_", gsub("ICR_([0-9]+).*","\\1",imp_whole$ID))
  # imp_whole$icr_name <- str_replace( imp_whole$icr_id, "_", " ")
  # Scan for previously published icrs
  high_conf_icrs <- imp_whole$icr_id[grep('#', imp_whole$ID)]

  # Get list of ICRs that have evidence of gametic origin for methylation
  # Includes "high confidence" (lit validated ICRs)
  imp_gamete <- tdhia::imprintome_icr_gametic_nearest_transcripts
  imp_gamete$icr_id <- paste0("ICR_", gsub("ICR_([0-9]+).*","\\1",imp_gamete$ID))
  # med_conf_icrs <- imp_gamete$icr_id #setdiff(imp_gamete$icr_id, high_conf_icrs)
  med_conf_icrs <- tdhia::imprintome_gametic_icrs$icr_id
  low_conf_icrs <- setdiff(imp_whole$icr_id, union(med_conf_icrs, high_conf_icrs))


  # Add zinc finger info
  df$is_icr_zinc <- df$icr_id %in% df_zinc_finger$icr_id[df_zinc_finger$near_zinf_finger]
  df$ZFP57  <- df$icr_id %in% df_zinc_finger$icr_id[df_zinc_finger$ZFP57]
  df$ZFP445 <- df$icr_id %in% df_zinc_finger$icr_id[df_zinc_finger$ZFP445]
  
  # Add icr confidence info
  df$icr_conf = 3*(df$icr_id %in% low_conf_icrs)
  df$icr_conf[df$icr_id %in% med_conf_icrs] = 2
  df$icr_conf[df$icr_id %in% high_conf_icrs] = 1

  # Add closest gene to each ICR
  # df_gene <- merge(x = df, y = imp_whole %>%
  #                       dplyr::select("icr_id", "Genomic.Coordinates", "Nearest.Transcript",
  #                                     "Distance.to.Nearest.Transcript"),
  #                     by = "icr_id", all.x = TRUE, all.y = FALSE, sort = FALSE)

  df_gene <- dplyr::left_join(x = df, y = imp_whole %>%
              dplyr::select("icr_id", "Genomic.Coordinates", "Nearest.Transcript",
                            "Distance.to.Nearest.Transcript"), by = "icr_id", keep =  NULL)
  

  return(df_gene)
}


#' #' add_metadata_from_cpg
#' #'
#' #' @param cpg_id character vector of cpg ids
#' #' @importFrom rlang .data
#' #' @export
#' #' @returns dataframe with ICR ids, ICR confidence, and zinc finger status
#' add_metadata_from_cpg <- function(cpg_id) {
#'   
#'   
#'   df = data.frame(cpg_id = cpg_id)
#'   df$icr_id <- sapply(df$cpg_id, function(x) 
#'     tdhia::mapping_cpg_icr_ids$ICR_id[x==tdhia::mapping_cpg_icr_ids$CpG_id][1])
#'   
#'   
#'   # Read list of ICRs that overlap with nearby zing finger
#'   df_zinc_finger <- tdhia::imprintome_icr_zinc_finger
#'   df_zinc_finger$icr_id <- as.numeric(gsub("ICR_([0-9]+).*","\\1",df_zinc_finger$icr))
#'   df_zinc_finger$near_zinf_finger <- rowSums(!df_zinc_finger[, 3:5] == "", na.rm = TRUE) > 0
#'   zinc_finger_icrs <- df_zinc_finger$icr_id[df_zinc_finger$near_zinf_finger]
#'   
#'   
#'   # Load annotated list of whole imprintome
#'   # ICR_IDs that end with "#" are high confidence
#'   imp_whole <- tdhia::imprintome_icr_nearest_transcripts
#'   imp_whole$icr_id <- as.numeric(gsub("ICR_([0-9]+).*","\\1",imp_whole$ID))
#'   imp_whole$icr_name <- paste0("ICR ", imp_whole$icr_id)
#'   # Scan for previously published icrs
#'   high_conf_icrs <- imp_whole$icr_id[grep('#', imp_whole$ID)]
#'   
#'   # Get list of ICRs that have evidence of gametic origin for methylation
#'   # Includes "high confidence" (lit validated ICRs)
#'   imp_gamete <- tdhia::imprintome_icr_gametic_nearest_transcripts
#'   imp_gamete$icr_id <- as.numeric(gsub("ICR_([0-9]+).*","\\1",imp_gamete$ID))
#'   med_conf_icrs <- imp_gamete$icr_id #setdiff(imp_gamete$icr_id, high_conf_icrs)
#'   low_conf_icrs <- setdiff(imp_whole$icr_id, union(med_conf_icrs, high_conf_icrs))
#'   
#'   df  = data.frame(cpg_id = cpg_id)
#'  
#'   df2 <- left_join(x = df, y =  distinct(tdhia::manifest_v1A2_design_scores %>%
#'                                            select(cpg_id, icr_id)),
#'                    by = join_by(cpg_id == cpg_id), keep = FALSE, multiple = "first",
#'                    unmatched = "drop", na_matches = "never")
#'   # Filter entries without a valid ICR
#'   df2 <- df2 %>% filter(grepl(pattern = "^ICR_.*$", x = icr_id))
#'   
#'   df2$icr_conf = 3*(df2$icr_id %in% paste0("ICR_", low_conf_icrs))
#'   df2$icr_conf[df2$icr_id %in% paste0("ICR_", med_conf_icrs)] = 2
#'   df2$icr_conf[df2$icr_id %in% paste0("ICR_", high_conf_icrs)] = 1
#'   
#'   
#'   df2$near_zinc_finger = df2$icr_id %in% zinc_finger_icrs
#'   
#'   return(df2)
#' }