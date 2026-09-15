

#' Aggregate Probe Methylation to CpG Sites
#'
#' Averages beta values from probes assigned to the same CpG, with optional
#' mapping filters, quantile normalization, sorting, and within-ICR smoothing.
#'
#' @param probe_beta List containing probe_beta_df (a numeric data frame with
#'   probes in rows and samples in columns), platform, and manifest. The
#'   manifest must contain Name and MAPINFO for genomic mapping checks.
#'   Probe p-values are not used by this function.
#' @param quantile_norm Logical; quantile-normalize sample columns of the
#'   aggregated CpG beta data using preprocessCore::normalize.quantiles().
#' @param discard_unmapped_cpgs Logical; remove CpGs whose first matching
#'   manifest entry has MAPINFO equal to zero or NA.
#' @param discard_non_icr_cpgs Logical; remove CpGs without a matching ICR in
#'   mapping_cpg_icr_ids.
#' @param smooth_adj_cpgs Logical; sort by ICR and genomic position, then
#'   apply a centered rolling mean of width three within each ICR, separately
#'   for each sample. Partial windows are used at region boundaries.
#' @param sort_cpgs Logical; sort by numeric ICR suffix and CpG_start.
#'   Smoothing also performs this sorting regardless of sort_cpgs.
#' @param db_flag Logical; save the initial environment to
#'   convert_probes_to_cpgs.RData in the working directory.
#' @param multicore Number of workers for smoothing. Defaults to the detected
#'   core count minus two, with a minimum of one. Unused without smoothing.
#'
#' @details
#' CpG IDs are obtained by removing an underscore followed by four terminal
#' characters from probe row names. Means ignore missing beta measurements;
#' an all-missing group produces NaN. Probe counts include all contributing
#' probe rows, rather than only nonmissing values. Quantile normalization,
#' when requested, precedes sorting and smoothing. Progress messages are
#' printed unconditionally.
#'
#' @return A named list with:
#'   - cpg_beta_df: data frame of aggregated beta values, CpG IDs in rows
#'     and sample IDs in columns.
#'   - platform and manifest: metadata copied from probe_beta.
#'   - cpg_n_probes: data frame containing n_probes for retained CpGs.
#'     Match by CpG row names: this table is captured before optional sorting
#'     and smoothing and is not reordered with cpg_beta_df.
#'   - input_args: argument values excluding probe_beta.
#' @seealso [filter_probes()], [convert_cpgs_to_icrs()]
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
convert_probes_to_cpgs <- function(probe_beta, quantile_norm = FALSE,
                                   discard_unmapped_cpgs = TRUE,
                                   discard_non_icr_cpgs = TRUE,
                                   smooth_adj_cpgs = FALSE,
                                   sort_cpgs = FALSE, db_flag = FALSE,
                                   multicore = max(c(parallel::detectCores()-2,1))) {
  # Get manifest data
  mft = probe_beta$manifest
  icr_mapping = tdhia::mapping_cpg_icr_ids

  if (db_flag) {save(list = ls(all.names = TRUE), file = "convert_probes_to_cpgs.RData")}
  # load(file = "convert_probes_to_cpgs.RData")


  # Match CpG site ID for each probe ID in probe_beta_df2
  # probe_beta_df2 <-
  #   probe_beta$probe_beta_df %>%
  #   tibble::rownames_to_column(var = "Probe_ID") %>%
  #   dplyr::left_join(y = dplyr::select(mft, c("Probe_ID", "Name")),
  #                    "Probe_ID",-"Probe_ID") %>%
  #   dplyr::select(-c("Probe_ID")) %>%
  #   dplyr::rename("CpG_ID" = "Name")

  probe_beta_df2 <- probe_beta$probe_beta_df
  probe_beta_df2$CpG_ID = stringr::str_replace(rownames(probe_beta_df2),"_.{4}$","")

  # CpG Beta Matrix: average beta values between probe_id(s) that belong same CpG site
  cpg_beta_df <-
    probe_beta_df2 %>%
    dplyr::group_by(.data$CpG_ID) %>%
    dplyr::summarize(dplyr::across(dplyr::where(base::is.numeric),
                                   function(x) mean(x, na.rm = TRUE)),
                     n_probes = dplyr::n())
  cat(sprintf("CpG filter: %.0f of probes in dataset mapped to %.0f unique CpG sites id data.\n",
              nrow(probe_beta$probe_beta_df), nrow(cpg_beta_df)))

  # Convert to dataframe to add rownames and drop first column (Cpg IDs)
  cpg_beta_df <- as.data.frame(cpg_beta_df)
  rownames(cpg_beta_df) <- cpg_beta_df$CpG_ID
  cpg_beta_df <- dplyr::select(cpg_beta_df,-.data$CpG_ID)


  # Discard CpG_ID(s) that do not map uniquely to a single genome location
  #   Defined in manifest file where MAPINFO == (0 or NA)
  if (discard_unmapped_cpgs) {
    f_unmmaped = function (x) base::is.na(x) | x==0;
    unmapped_cpg_ids <- sapply(rownames(cpg_beta_df),
                               function(x) f_unmmaped(mft$MAPINFO[which(x==mft$Name)[1]]))

    # Remove unmapped cpgs
    cpg_beta_df <- cpg_beta_df[!unmapped_cpg_ids,]

    cat(sprintf("CpG filter: discarded %.0f%% of CpG sites (%i/ %i) because they
      do not map uniquely to the genome. %i CpG sites remain.\n",
                100*sum(unmapped_cpg_ids)/length(unmapped_cpg_ids), sum(unmapped_cpg_ids),
                length(unmapped_cpg_ids), nrow(cpg_beta_df)))

  }


  if (discard_non_icr_cpgs){
    discard_non_icr <- cpg_beta_df %>%
      tibble::rownames_to_column("CpG_id") %>%
      dplyr::left_join(y = dplyr::select(icr_mapping, c("CpG_id", "ICR_id")),
                       by = "CpG_id",unmatched = "drop",
                       multiple = "first") %>%
      dplyr::pull("ICR_id") %>% base::is.na()
    # Remove non-icr cpgs
    cpg_beta_df <- cpg_beta_df[!discard_non_icr,]

    cat(sprintf("CpG filter: discarded %.0f%% of CpG sites (%i/ %i) because they
    do not map to an ICR. %i CpG sites remain.\n",
                100*sum(discard_non_icr)/nrow(cpg_beta_df), sum(discard_non_icr),
                nrow(cpg_beta_df), nrow(cpg_beta_df)))
  }


  # Extract n_probes from dataframe (to be stored separately)
  n_probes <- dplyr::select(cpg_beta_df, n_probes)
  cpg_beta_df <- cpg_beta_df %>% dplyr::select(-n_probes)

  if (quantile_norm) {
    # Quantile Normalization by column
    cat("Performing quantile normalization...\n")
    norm_mat <- preprocessCore::normalize.quantiles(as.matrix(cpg_beta_df))
    dimnames(norm_mat) <- list(rownames(cpg_beta_df), colnames(cpg_beta_df))
    cpg_beta_df <- norm_mat
  }


  # Sorting function used below
  cpg_sort_fun <- function(cpg_beta_df) {# Order cpg_beta
    cpg_mappings <- cpg_beta_df %>%
      tibble::rownames_to_column("CpG_id") %>%
      dplyr::left_join(y = dplyr::select(icr_mapping, c("CpG_id", "ICR_id",
                                                        "CpG_start")),
                       by = "CpG_id",unmatched = "drop",
                       multiple = "first")
    # Sort by ICR id, then cpg position
    cpg_mappings <-
      cpg_mappings[order(as.numeric(stringr::str_replace(cpg_mappings$ICR_id, "ICR_","")),
                         cpg_mappings$CpG_start,decreasing = FALSE),]
    rownames(cpg_mappings) <- cpg_mappings$CpG_id

    return(cpg_mappings)
  }


  # Sort cpgs by ICR id and genomic position
  if (sort_cpgs) {
    cat("Cpg Filter: Sorting CpGs by ICR and genomic index.\n")
    cpg_beta_df <-  dplyr::select(cpg_sort_fun(cpg_beta_df),-c("CpG_id", "ICR_id",
                                                               "CpG_start"))
  }


  # Perform rolling average on beta values with adjacent cpgs
  if (smooth_adj_cpgs) {
    cat("CpG Filter: Sorting and smoothing adjacent cpg beta values.\n")

    # Sort cpgs
    cpg_mappings <- cpg_sort_fun(cpg_beta_df)

    # Get unique ICR IDs
    unq_icr_ids <- unique(cpg_mappings$ICR_id)

    # Parallel processing
    cl <- parallel::makeCluster(multicore)
    doParallel::registerDoParallel(cl)


    # Grab all cpgs for each ICR, do column wise smoothing
    n=NULL # Prevent a devtools::check() warning
    sm_cpg_df <- foreach::foreach(
      n = seq_along(unq_icr_ids), .combine = 'rbind',
      .packages = c("dplyr","zoo")) %dopar% {
        ix <- unq_icr_ids[n]==cpg_mappings$ICR_id
        cpg_sub <- cpg_mappings[unq_icr_ids[n]==cpg_mappings$ICR_id,] %>%
          dplyr::select(-c("CpG_id", "ICR_id", "CpG_start"))
        rownames(cpg_sub)
        
        sm_cpg_sub <-
          zoo::rollapply(cpg_sub, width=3, FUN = function(x) mean(x, na.rm = TRUE),
                         by = 1, by.column = TRUE, fill = NA, align ="center",
                         partial = TRUE)
        rownames(sm_cpg_sub)  <- cpg_mappings$CpG_id[unq_icr_ids[n] == cpg_mappings$ICR_id]
        sm_cpg_sub <- as.data.frame(sm_cpg_sub)
        
        sm_cpg_sub
      }
    #stop cluster
    parallel::stopCluster(cl)

    # Overwrite
    cpg_beta_df <- sm_cpg_df

  }

  # Export data
  cpg_beta <- list(
    cpg_beta_df = as.data.frame(cpg_beta_df), platform = probe_beta$platform,
    manifest = probe_beta$manifest, cpg_n_probes = n_probes,
    input_args = mget(setdiff(names(formals()), c("probe_beta")), envir = environment()))

  return(cpg_beta)
}
