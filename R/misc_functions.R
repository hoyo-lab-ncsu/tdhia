



#' match_cpg_id_to_icr_id
#'
#' @description looks up corresponding ICR ids given a vector of cpg_ids
#'
#' @param cpg_ids vector of cpg_ids
#' @param AsNumeric flag, when true, returns in numeric form
#' @return vector of ICR ids
#' @importFrom rlang .data
#' @export
match_cpg_id_to_icr_id <- function(cpg_ids, AsNumeric = FALSE) {

  df_mapping <- tdhia::mapping_cpg_icr_ids

  df_icr <- dplyr::left_join(x = data.frame(CpG_id = cpg_ids), y = dplyr::distinct(dplyr::select(
    df_mapping, c(.data$CpG_id, .data$ICR_id))), keep = NULL, by = "CpG_id")

  if (AsNumeric) {
    out <- as.numeric(stringr::str_replace(df_icr$ICR_id, "ICR_", ""))
  } else {
    out <- df_icr$ICR_id
  }

  return(out)
}
