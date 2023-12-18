
#' load_idata_to_probes
#'
#' Loads a series flourescence measurements from a custom infinium methylation array
#' and converts into methylation beta values.
#'
#' beta = M/(M+U)   (for probes with significant detection p-values)
#'
#' @param idat_dir_paths the full filesystem path to directory containing the idat
#' files to be processed, formatted as a string or a vector of strings.
#' @param platform string for the platform that the array belongs to, as specifed
#' within the SeSame package with openSesame().
#' @param mft manifest file for particular imprintome array used in data collection.
#' A dataframe with metadata mapping probe_ids to CpG sites and genome locations.
#' See ?manifest_v1A2 for more info.
#' @param multicore boolean flag, when true the max number of cores minus 1 is
#' used for the current system.
#' @param idat_name_remapping a dataframe containing two columns used to convert the
#' basename of the IDAT files (used to label columns of output dataframe), into
#' a different set of identifiers/mappings.
#'  - idat_basename: column of base filenames of IDAT (no repeated values).
#'  - id: column of strings of new mapping (no repeated values).
#' @param idat_discard_unmapped a boolean flag, when set to TRUE discards any idat base
#'  filenames that are not included in the idat_name_remapping argument (note:
#'  idat_name_remapping must be specified).
#' @param quantile_norm a boolean flag, when set to TRUE, applies quantile
#' normalization between the columns in the probe_beta dataframe (default = FALSE).
#' @param mask a boolean flag specified in opensesame() in the Sesame package,
#' when TRUE excludes some probes due to issues of inter-dependence of measurements
#' (default = FALSE).
#'
#' @returns a named list object with the following fields:
#'  - probe_beta_df: dataframe of beta values, probe_id x sample_id
#'  - probe_pval_df: dataframe of signal p-values, probe_id x sample_id
#'  - platform: string that describes platform for methylation array
#'  - manifest: dataframe of the manifest file used for SeSame processing
#'
load_idata_to_probes <- function(idat_dir_paths, platform = "TruDx_imprintome", mft = NULL, multicore = TRUE, idat_name_remapping = NULL, idat_discard_unmapped = TRUE,
    quantile_norm = FALSE, mask = FALSE) {
    # Load manifest file if platform is ture diagnostic imprintome array
    if (is.null(mft) && platform == "TruDx_imprintome") {
        mft = tdhia::manifest_v1A2
    }

    # save(list = ls(all.names = TRUE), file = 'load_idata_debug.RData') load(file = 'load_idata_debug.RData')



    # Get list of IDAT files in target path
    all_idat_names <- dir(path = idat_dir_paths, pattern = "*.idat", ignore.case = TRUE)
    # Extract the basenames of IDAT files Removes _Grn.idat and _Red.idat from file names (case insensitive)
    all_idat_basenames <- stringr::str_replace(all_idat_names, pattern = stringr::regex("_(Grn|Red).idat$", ignore_case = TRUE), replacement = "")


    # Debugging: verify that two IDAT files are associated with each basename: One file for red and one for green channel fluorescence
    basename_tbl <- as.data.frame(table(all_idat_basenames))
    if (!all(basename_tbl$Freq == 2)) {
        warning(paste("IDAT_PAIR: Not all IDAT file basenames have 2 IDAT files associated with them.", "Files are missing! See basename_tbl variable"))
    }

    # Get unique list of basenames (removes repeats)
    idat_basenames <- unique(all_idat_basenames)

    # Discard IDAT files that are not included in idat_name_remapping file
    if (!is.null(idat_name_remapping) && idat_discard_unmapped) {
        idat_basenames <- idat_basenames[is.element(idat_basenames, idat_name_remapping$idat_basename)]
    }


    # Caching SeSame Data Files (required for running SeSame firs time...)
    cat("Checking that Sesame idatSignature file is locally cached...")
    sesameData::sesameDataCache()
    cat("Done.\n")

    # Set Multicore parameter depending on operating system
    default_ncores <- max(c(parallel::detectCores() - 1, 1))
    if (.Platform$OS.type == "windows" && multicore) {
        multicore_arg <- BiocParallel::SnowParam(default_ncores)
    } else if (multicore) {
        multicore_arg <- BiocParallel::MulticoreParam(default_ncores)
    } else {
        multicore_arg = NULL
    }

    # Load each of the IDAT file pairs, process with standard SeSame pipeline
    cat("Processing IDAT Files...")
    if (!is.null(multicore_arg)) {
        cat(sprintf(" (Using %i/%i cores.)...", default_ncores, parallel::detectCores()))
    }
    probe_beta_matrix <- sesame::openSesame(paste0(idat_dir_paths, "/", idat_basenames), platform = platform, manifest = mft, mask = mask, BPPARAM = multicore_arg,
        fun = sesame::getBetas)


    probe_pval_matrix <- sesame::openSesame(paste0(idat_dir_paths, "/", idat_basenames), platform = platform, manifest = mft, BPPARAM = multicore_arg,
        fun = sesame::pOOBAH, return.pval = TRUE)

    # Debugging: saving output from sesame pipeline Todo: pipeline is currently run twice to get beta and p-values, recode to run once
    # save(list = ls(all.names = TRUE), file = 'load_idata_debug.RData') load(file = 'load_idata_debug.RData')


    # Convert matrix to dataframe for easier manipulation
    probe_beta_df <- as.data.frame(probe_beta_matrix)
    probe_pval_df <- as.data.frame(probe_pval_matrix)
    cat("Done.\n")

    # If a column mapping is specified, rename column names in probe_beta_df Discard unmapped columns if user specifies it
    if (!is.null(idat_name_remapping)) {

        # If idat_discard_unmapped is TRUE, any columns in [probe_beta_df] not included in [column mapping] are discarded. This is useful if
        # the IDAT folder has data from multiple studies
        if (idat_discard_unmapped) {
            probe_beta_df <- probe_beta_df[, is.element(colnames(probe_beta_df), idat_name_remapping$idat_basename)]
            probe_pval_df <- probe_pval_df[, is.element(colnames(probe_pval_df), idat_name_remapping$idat_basename)]
        }
        # Rename by matching new id to each idat basename using a left_join()
        remapped_colnames <- dplyr::left_join(df_remap <- data.frame(idat_basename = colnames(probe_beta_df)), dplyr::select(idat_name_remapping,
            c("idat_basename", "id")), by = "idat_basename", na_matches = "never", unmatched = "error", relationship = "one-to-one")
        # Rename columns of probe_beta_df and probe_pval_df
        colnames(probe_beta_df) <- remapped_colnames$id
        colnames(probe_pval_df) <- remapped_colnames$id
    }


    # Apply quantile normalization between all columns if specified by user.
    if (quantile_norm) {
        norm_probe_df <- as.data.frame(preprocessCore::normalize.quantiles(as.matrix(probe_beta_df)), row.names = rownames(probe_beta_df))
        colnames(norm_probe_df) <- colnames(probe_beta_df)
        probe_beta_df <- norm_probe_df
    }

    # # Add probe_id as a first column and remove row names probe_beta_df<-cbind(data.frame(Probe_ID =
    # rownames(probe_beta_df)),probe_beta_df) rownames(probe_beta_df) <- NULL probe_pval_df<-cbind(data.frame(Probe_ID =
    # rownames(probe_pval_df)),probe_pval_df) rownames(probe_pval_df) <- NULL



    probe_beta <- list(probe_beta_df = probe_beta_df, probe_pval_df = probe_pval_df, platform = platform, manifest = mft)
    return(probe_beta)

}


