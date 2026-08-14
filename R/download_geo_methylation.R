

#' Download methylation-array IDAT files and GEO metadata
#'
#' Downloads IDAT files associated with a GEO Series accession and creates
#' sample metadata CSV files from the GEO Series Matrix records. Files are
#' stored in a persistent user-data directory associated with tdhia.
#'
#' Existing files are not downloaded again unless `overwrite = TRUE`.
#'
#' @param geo_accession Character scalar. A GEO Series accession such as
#'   `"GSE268075"`.
#' @param data_dir Character scalar. Root directory in which GEO data are
#'   stored. Defaults to the package-specific persistent user-data directory
#'   returned by `tools::R_user_dir("tdhia", "data")`.
#' @param overwrite Logical. If `TRUE`, download and recreate existing files.
#' @param extract_archives Logical. If `TRUE`, extract `.tar`, `.tar.gz`, and
#'   `.tgz` archives that may contain IDAT files.
#' @param decompress_idats Logical. If `TRUE`, decompress `.idat.gz` files
#'   in the IDAT directory and remove the compressed copies after successful
#'   decompression.
#' @param remove_archives Logical. If `TRUE`, remove downloaded archives after
#'   successful extraction. Defaults to `FALSE`.
#' @param quiet Logical. Suppress download progress messages where possible.
#'
#' @return Invisibly returns a list containing:
#' \describe{
#'   \item{geo_accession}{Normalized GEO accession.}
#'   \item{directory}{Dataset directory.}
#'   \item{idat_files}{Paths to discovered, decompressed IDAT files.}
#'   \item{metadata_files}{Paths to generated metadata CSV files.}
#'   \item{supplementary_files}{Paths to downloaded GEO supplementary files.}
#' }
#'
#' @details
#' GEO commonly distributes methylation-array IDAT files either as individual
#' `.idat` or `.idat.gz` files or inside archives such as a `RAW.tar` file.
#' This function handles both arrangements.
#'
#' Metadata are obtained from the GEO Series Matrix using
#' `GEOquery::getGEO()` and written as one CSV file per platform.
#'
#' @seealso [tdhia_geo_data_dir()], [tdhia_clear_geo_data()]
#'
#' @examples result <- tdhia::download_geo_methylation("GSE268075")
#' @author ChatGPT OpenAI GPT-5.6, tested by Bruce Corliss
#' @export
download_geo_methylation <- function(
    geo_accession,
    data_dir = tools::R_user_dir("tdhia", which = "data"),
    overwrite = FALSE,
    extract_archives = TRUE,
    decompress_idats = TRUE,
    remove_archives = FALSE,
    quiet = FALSE
) {
    if (!requireNamespace("GEOquery", quietly = TRUE)) {
        stop(
            "Package 'GEOquery' is required. Install it with:\n",
            "BiocManager::install('GEOquery')",
            call. = FALSE
        )
    }

    if (!requireNamespace("Biobase", quietly = TRUE)) {
        stop(
            "Package 'Biobase' is required. Install it with:\n",
            "BiocManager::install('Biobase')",
            call. = FALSE
        )
    }

    geo_accession <- validate_geo_series_accession(geo_accession)

    checkmate_flag <- function(x, name) {
        if (!is.logical(x) || length(x) != 1L || is.na(x)) {
            stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
        }
    }

    checkmate_flag(overwrite, "overwrite")
    checkmate_flag(extract_archives, "extract_archives")
    checkmate_flag(decompress_idats, "decompress_idats")
    checkmate_flag(remove_archives, "remove_archives")
    checkmate_flag(quiet, "quiet")

    if (!is.character(data_dir) ||
        length(data_dir) != 1L ||
        is.na(data_dir) ||
        !nzchar(data_dir)) {
        stop("`data_dir` must be a non-empty character scalar.", call. = FALSE)
    }

    dataset_dir <- file.path(data_dir, "geo", geo_accession)
    supplementary_dir <- file.path(dataset_dir, "supplementary")
    idat_dir <- file.path(dataset_dir, "idat")
    metadata_dir <- file.path(dataset_dir, "metadata")
    geoquery_cache_dir <- file.path(dataset_dir, "geoquery_cache")

    dirs <- c(
        dataset_dir,
        supplementary_dir,
        idat_dir,
        metadata_dir,
        geoquery_cache_dir
    )

    for (directory in dirs) {
        if (!dir.exists(directory) &&
            !dir.create(directory, recursive = TRUE, showWarnings = FALSE)) {
            stop("Could not create directory: ", directory, call. = FALSE)
        }
    }

    message_if_verbose(
        quiet,
        "Using GEO data directory: ",
        normalizePath(dataset_dir, winslash = "/", mustWork = FALSE)
    )

    supplementary_files <- download_geo_idat_supplementary_files(
        geo_accession = geo_accession,
        destination = supplementary_dir,
        overwrite = overwrite,
        quiet = quiet
    )

    if (extract_archives) {
        extract_geo_idat_archives(
            files = supplementary_files,
            destination = idat_dir,
            overwrite = overwrite,
            remove_archives = remove_archives,
            quiet = quiet
        )
    }

    copy_loose_idat_files(
        source_directory = supplementary_dir,
        destination = idat_dir,
        overwrite = overwrite,
        quiet = quiet
    )

    if (decompress_idats) {
        decompress_geo_idat_files(
            directory = idat_dir,
            overwrite = overwrite,
            remove_compressed = TRUE,
            quiet = quiet
        )
    }

    metadata_files <- download_geo_sample_metadata(
        geo_accession = geo_accession,
        metadata_dir = metadata_dir,
        geoquery_cache_dir = geoquery_cache_dir,
        overwrite = overwrite,
        quiet = quiet
    )

    idat_pattern <- if (decompress_idats) {
        "\\.idat$"
    } else {
        "\\.idat(?:\\.gz)?$"
    }

    idat_files <- list.files(
        idat_dir,
        pattern = idat_pattern,
        recursive = TRUE,
        full.names = TRUE,
        ignore.case = TRUE
    )

    idat_files <- sort(unique(normalizePath(
        idat_files,
        winslash = "/",
        mustWork = FALSE
    )))

    if (!length(idat_files)) {
        warning(
            "No IDAT files were found for ", geo_accession, ". ",
            "The GEO entry may not provide raw IDAT files, or its files may ",
            "use an unexpected archive format.",
            call. = FALSE
        )
    }

    message_if_verbose(
        quiet,
        "Found ",
        length(idat_files),
        " IDAT file",
        if (length(idat_files) == 1L) "" else "s",
        "."
    )

    invisible(list(
        geo_accession = geo_accession,
        directory = normalizePath(
            dataset_dir,
            winslash = "/",
            mustWork = FALSE
        ),
        idat_files = idat_files,
        metadata_files = metadata_files,
        supplementary_files = supplementary_files
    ))
}


validate_geo_series_accession <- function(geo_accession) {
    if (!is.character(geo_accession) ||
        length(geo_accession) != 1L ||
        is.na(geo_accession) ||
        !nzchar(geo_accession)) {
        stop(
            "`geo_accession` must be a non-empty character scalar.",
            call. = FALSE
        )
    }

    geo_accession <- toupper(trimws(geo_accession))

    if (!grepl("^GSE[0-9]+$", geo_accession)) {
        stop(
            "`geo_accession` must be a GEO Series accession such as ",
            "'GSE123456'.",
            call. = FALSE
        )
    }

    geo_accession
}


message_if_verbose <- function(quiet, ...) {
    if (!quiet) {
        message(...)
    }

    invisible(NULL)
}


download_geo_idat_supplementary_files <- function(
    geo_accession,
    destination,
    overwrite,
    quiet
) {
    message_if_verbose(
        quiet,
        "Retrieving supplementary-file listing for ",
        geo_accession,
        "..."
    )

    listing <- GEOquery::getGEOSuppFiles(
        GEO = geo_accession,
        fetch_files = FALSE
    )

    if (!is.data.frame(listing) || !nrow(listing)) {
        warning(
            "No supplementary files were listed for ",
            geo_accession,
            ".",
            call. = FALSE
        )
        return(character())
    }

    supplementary_table <- standardize_geo_supplementary_listing(listing)

    keep <- grepl(
        paste0(
            "(?i)",
            "\\.idat(?:\\.gz)?$",
            "|raw.*\\.(?:tar|tar\\.gz|tgz)$",
            "|\\.(?:tar|tar\\.gz|tgz)$"
        ),
        supplementary_table$filename,
        perl = TRUE
    )

    supplementary_table <- supplementary_table[keep, , drop = FALSE]

    if (!nrow(supplementary_table)) {
        warning(
            "No IDAT files or candidate IDAT archives were found in the ",
            "supplementary-file listing for ",
            geo_accession,
            ".",
            call. = FALSE
        )
        return(character())
    }

    output_files <- file.path(
        destination,
        basename(supplementary_table$filename)
    )

    for (i in seq_len(nrow(supplementary_table))) {
        output_file <- output_files[[i]]
        url <- supplementary_table$url[[i]]

        if (file.exists(output_file) && !overwrite) {
            message_if_verbose(
                quiet,
                "Already downloaded: ",
                basename(output_file)
            )
            next
        }

        message_if_verbose(
            quiet,
            "Downloading: ",
            basename(output_file)
        )

        temporary_file <- paste0(output_file, ".partial")

        if (file.exists(temporary_file)) {
            unlink(temporary_file, force = TRUE)
        }

        status <- tryCatch(
            utils::download.file(
                url = url,
                destfile = temporary_file,
                mode = "wb",
                quiet = quiet
            ),
            error = function(e) {
                unlink(temporary_file, force = TRUE)
                stop(
                    "Failed to download ",
                    basename(output_file),
                    ":\n",
                    conditionMessage(e),
                    call. = FALSE
                )
            }
        )

        if (!identical(status, 0L) || !file.exists(temporary_file)) {
            unlink(temporary_file, force = TRUE)
            stop(
                "Download did not complete successfully: ",
                basename(output_file),
                call. = FALSE
            )
        }

        if (file.exists(output_file)) {
            unlink(output_file, force = TRUE)
        }

        if (!file.rename(temporary_file, output_file)) {
            copied <- file.copy(
                temporary_file,
                output_file,
                overwrite = TRUE
            )

            unlink(temporary_file, force = TRUE)

            if (!copied) {
                stop(
                    "Could not move downloaded file into place: ",
                    output_file,
                    call. = FALSE
                )
            }
        }
    }

    output_files[file.exists(output_files)]
}


standardize_geo_supplementary_listing <- function(listing) {
    column_names <- names(listing)
    lower_names <- tolower(column_names)

    url_column <- which(
        grepl("url|ftp|download", lower_names)
    )

    filename_column <- which(
        grepl("file|name", lower_names)
    )

    if (!length(url_column)) {
        row_values <- rownames(listing)

        if (!is.null(row_values) &&
            all(grepl("^(?:https?|ftp)://", row_values))) {
            urls <- row_values
        } else {
            stop(
                "Could not identify the URL column returned by ",
                "GEOquery::getGEOSuppFiles().",
                call. = FALSE
            )
        }
    } else {
        urls <- as.character(listing[[url_column[[1L]]]])
    }

    if (length(filename_column)) {
        filenames <- as.character(
            listing[[filename_column[[1L]]]]
        )
    } else {
        filenames <- basename(sub("\\?.*$", "", urls))
    }

    missing_filename <- is.na(filenames) | !nzchar(filenames)

    if (any(missing_filename)) {
        filenames[missing_filename] <- basename(
            sub("\\?.*$", "", urls[missing_filename])
        )
    }

    data.frame(
        url = urls,
        filename = filenames,
        stringsAsFactors = FALSE
    )
}


extract_geo_idat_archives <- function(
    files,
    destination,
    overwrite,
    remove_archives,
    quiet
) {
    archives <- files[
        grepl(
            "\\.(?:tar|tar\\.gz|tgz)$",
            files,
            ignore.case = TRUE
        )
    ]

    archives <- archives[file.exists(archives)]

    if (!length(archives)) {
        return(invisible(character()))
    }

    extracted_files <- character()

    for (archive in archives) {
        message_if_verbose(
            quiet,
            "Inspecting archive: ",
            basename(archive)
        )

        archive_contents <- tryCatch(
            utils::untar(archive, list = TRUE),
            error = function(e) {
                warning(
                    "Could not inspect archive ",
                    basename(archive),
                    ": ",
                    conditionMessage(e),
                    call. = FALSE
                )
                character()
            }
        )

        idat_members <- archive_contents[
            grepl(
                "\\.idat(?:\\.gz)?$",
                archive_contents,
                ignore.case = TRUE
            )
        ]

        if (!length(idat_members)) {
            message_if_verbose(
                quiet,
                "No IDAT files found in ",
                basename(archive),
                "."
            )
            next
        }

        expected_files <- file.path(
            destination,
            basename(idat_members)
        )

        all_present <- all(file.exists(expected_files))

        if (all_present && !overwrite) {
            message_if_verbose(
                quiet,
                "IDAT files from ",
                basename(archive),
                " have already been extracted."
            )

            extracted_files <- c(extracted_files, expected_files)
            next
        }

        extraction_directory <- tempfile(
            pattern = paste0(
                "tdhia_",
                tools::file_path_sans_ext(basename(archive)),
                "_"
            )
        )

        dir.create(
            extraction_directory,
            recursive = TRUE,
            showWarnings = FALSE
        )

        extraction_ok <- FALSE

        tryCatch(
            {
                message_if_verbose(
                    quiet,
                    "Extracting ",
                    length(idat_members),
                    " IDAT file",
                    if (length(idat_members) == 1L) "" else "s",
                    " from ",
                    basename(archive),
                    "..."
                )

                utils::untar(
                    tarfile = archive,
                    files = idat_members,
                    exdir = extraction_directory
                )

                extraction_ok <- TRUE
            },
            error = function(e) {
                stop(
                    "Failed to extract ",
                    basename(archive),
                    ":\n",
                    conditionMessage(e),
                    call. = FALSE
                )
            }
        )

        if (!extraction_ok) {
            unlink(extraction_directory, recursive = TRUE, force = TRUE)
            next
        }

        extracted <- list.files(
            extraction_directory,
            pattern = "\\.idat(?:\\.gz)?$",
            recursive = TRUE,
            full.names = TRUE,
            ignore.case = TRUE
        )

        for (source_file in extracted) {
            destination_file <- file.path(
                destination,
                basename(source_file)
            )

            if (file.exists(destination_file) && !overwrite) {
                next
            }

            copied <- file.copy(
                from = source_file,
                to = destination_file,
                overwrite = overwrite
            )

            if (!copied) {
                unlink(extraction_directory, recursive = TRUE, force = TRUE)
                stop(
                    "Could not copy extracted IDAT file to: ",
                    destination_file,
                    call. = FALSE
                )
            }
        }

        extracted_files <- c(
            extracted_files,
            file.path(destination, basename(extracted))
        )

        unlink(
            extraction_directory,
            recursive = TRUE,
            force = TRUE
        )

        if (remove_archives) {
            unlink(archive, force = TRUE)
        }
    }

    invisible(unique(extracted_files[file.exists(extracted_files)]))
}


copy_loose_idat_files <- function(
    source_directory,
    destination,
    overwrite,
    quiet
) {
    loose_idats <- list.files(
        source_directory,
        pattern = "\\.idat(?:\\.gz)?$",
        recursive = TRUE,
        full.names = TRUE,
        ignore.case = TRUE
    )

    if (!length(loose_idats)) {
        return(invisible(character()))
    }

    copied_files <- character()

    for (source_file in loose_idats) {
        destination_file <- file.path(
            destination,
            basename(source_file)
        )

        if (file.exists(destination_file) && !overwrite) {
            message_if_verbose(
                quiet,
                "Already available: ",
                basename(destination_file)
            )

            copied_files <- c(copied_files, destination_file)
            next
        }

        copied <- file.copy(
            from = source_file,
            to = destination_file,
            overwrite = overwrite
        )

        if (!copied) {
            stop(
                "Could not copy IDAT file to: ",
                destination_file,
                call. = FALSE
            )
        }

        copied_files <- c(copied_files, destination_file)
    }

    invisible(copied_files)
}


decompress_geo_idat_files <- function(
    directory,
    overwrite = FALSE,
    remove_compressed = TRUE,
    quiet = FALSE
) {
    compressed_files <- list.files(
        path = directory,
        pattern = "\\.idat\\.gz$",
        recursive = TRUE,
        full.names = TRUE,
        ignore.case = TRUE
    )

    if (!length(compressed_files)) {
        message_if_verbose(
            quiet,
            "No compressed IDAT files need decompression."
        )

        return(invisible(character()))
    }

    output_files <- character(length(compressed_files))

    for (i in seq_along(compressed_files)) {
        compressed_file <- compressed_files[[i]]

        output_file <- sub(
            "\\.gz$",
            "",
            compressed_file,
            ignore.case = TRUE
        )

        output_files[[i]] <- output_file

        if (file.exists(output_file) && !overwrite) {
            message_if_verbose(
                quiet,
                "Uncompressed IDAT already exists: ",
                basename(output_file)
            )

            if (remove_compressed) {
                removed <- unlink(
                    compressed_file,
                    force = TRUE
                )

                if (!identical(removed, 0L)) {
                    warning(
                        "Could not remove compressed IDAT file: ",
                        compressed_file,
                        call. = FALSE
                    )
                }
            }

            next
        }

        message_if_verbose(
            quiet,
            "Decompressing: ",
            basename(compressed_file)
        )

        temporary_file <- paste0(output_file, ".partial")

        if (file.exists(temporary_file)) {
            unlink(temporary_file, force = TRUE)
        }

        input_connection <- gzfile(
            compressed_file,
            open = "rb"
        )

        output_connection <- file(
            temporary_file,
            open = "wb"
        )

        success <- FALSE

        tryCatch(
            {
                repeat {
                    buffer <- readBin(
                        input_connection,
                        what = "raw",
                        n = 1024L * 1024L
                    )

                    if (!length(buffer)) {
                        break
                    }

                    writeBin(
                        buffer,
                        output_connection
                    )
                }

                success <- TRUE
            },
            error = function(e) {
                stop(
                    "Failed to decompress ",
                    basename(compressed_file),
                    ":\n",
                    conditionMessage(e),
                    call. = FALSE
                )
            },
            finally = {
                close(input_connection)
                close(output_connection)

                if (!success && file.exists(temporary_file)) {
                    unlink(temporary_file, force = TRUE)
                }
            }
        )

        if (!file.exists(temporary_file) ||
            is.na(file.info(temporary_file)$size) ||
            file.info(temporary_file)$size == 0) {
            unlink(temporary_file, force = TRUE)

            stop(
                "Decompression produced an empty or missing file for: ",
                compressed_file,
                call. = FALSE
            )
        }

        if (file.exists(output_file)) {
            unlink(output_file, force = TRUE)
        }

        moved <- file.rename(
            temporary_file,
            output_file
        )

        if (!moved) {
            moved <- file.copy(
                from = temporary_file,
                to = output_file,
                overwrite = TRUE
            )

            if (moved) {
                unlink(temporary_file, force = TRUE)
            }
        }

        if (!moved || !file.exists(output_file)) {
            stop(
                "Could not move the decompressed IDAT file into place: ",
                output_file,
                call. = FALSE
            )
        }

        if (remove_compressed) {
            removed <- unlink(
                compressed_file,
                force = TRUE
            )

            if (!identical(removed, 0L)) {
                warning(
                    "The IDAT file was decompressed, but the compressed ",
                    "file could not be removed: ",
                    compressed_file,
                    call. = FALSE
                )
            }
        }
    }

    invisible(
        normalizePath(
            output_files[file.exists(output_files)],
            winslash = "/",
            mustWork = FALSE
        )
    )
}


download_geo_sample_metadata <- function(
    geo_accession,
    metadata_dir,
    geoquery_cache_dir,
    overwrite,
    quiet
) {
    existing_metadata <- list.files(
        metadata_dir,
        pattern = paste0(
            "^",
            geo_accession,
            "_metadata(?:_.+)?\\.csv$"
        ),
        full.names = TRUE,
        ignore.case = TRUE
    )

    if (length(existing_metadata) && !overwrite) {
        message_if_verbose(
            quiet,
            "Metadata already available for ",
            geo_accession,
            "."
        )

        return(sort(existing_metadata))
    }

    message_if_verbose(
        quiet,
        "Downloading GEO Series Matrix metadata for ",
        geo_accession,
        "..."
    )

    geo_objects <- GEOquery::getGEO(
        GEO = geo_accession,
        destdir = geoquery_cache_dir,
        GSEMatrix = TRUE,
        AnnotGPL = FALSE,
        getGPL = FALSE
    )

    if (!is.list(geo_objects)) {
        geo_objects <- list(geo_objects)
    }

    if (!length(geo_objects)) {
        stop(
            "No GEO Series Matrix records were returned for ",
            geo_accession,
            ".",
            call. = FALSE
        )
    }

    metadata_files <- character(length(geo_objects))

    for (i in seq_along(geo_objects)) {
        geo_object <- geo_objects[[i]]

        metadata <- Biobase::pData(geo_object)

        if (!is.data.frame(metadata)) {
            metadata <- as.data.frame(
                metadata,
                stringsAsFactors = FALSE
            )
        }

        metadata$geo_accession <- rownames(metadata)

        metadata <- metadata[
            ,
            c(
                "geo_accession",
                setdiff(names(metadata), "geo_accession")
            ),
            drop = FALSE
        ]

        platform <- tryCatch(
            Biobase::annotation(geo_object),
            error = function(e) ""
        )

        if (!length(platform) ||
            is.na(platform) ||
            !nzchar(platform)) {
            platform <- paste0("matrix", i)
        }

        platform <- gsub(
            "[^A-Za-z0-9._-]+",
            "_",
            platform
        )

        if (length(geo_objects) == 1L) {
            filename <- paste0(
                geo_accession,
                "_metadata.csv"
            )
        } else {
            filename <- paste0(
                geo_accession,
                "_metadata_",
                platform,
                ".csv"
            )
        }

        output_file <- file.path(metadata_dir, filename)

        utils::write.csv(
            metadata,
            file = output_file,
            row.names = FALSE,
            na = "",
            fileEncoding = "UTF-8"
        )

        metadata_files[[i]] <- output_file

        message_if_verbose(
            quiet,
            "Saved metadata: ",
            basename(output_file)
        )
    }

    sort(metadata_files)
}


#' Locate the persistent tdhia GEO data directory
#'
#' @param geo_accession Optional GEO Series accession. When supplied, returns
#'   the directory for that dataset.
#'
#' @return A character scalar containing the requested directory.
#'
#' @export
tdhia_geo_data_dir <- function(geo_accession = NULL) {
    directory <- file.path(
        tools::R_user_dir("tdhia", which = "data"),
        "geo"
    )

    if (!is.null(geo_accession)) {
        geo_accession <- validate_geo_series_accession(geo_accession)
        directory <- file.path(directory, geo_accession)
    }

    normalizePath(
        directory,
        winslash = "/",
        mustWork = FALSE
    )
}


#' Remove downloaded tdhia GEO data
#'
#' @param geo_accession Optional GEO Series accession. If supplied, only data
#'   for that accession are removed. If `NULL`, all GEO data downloaded by
#'   tdhia are removed.
#' @param confirm Logical. If `TRUE`, require interactive confirmation before
#'   deleting all downloaded GEO datasets.
#'
#' @return Invisibly returns `TRUE` if the directory no longer exists.
#'
#' @export
tdhia_clear_geo_data <- function(
    geo_accession = NULL,
    confirm = interactive()
) {
    if (!is.logical(confirm) ||
        length(confirm) != 1L ||
        is.na(confirm)) {
        stop("`confirm` must be TRUE or FALSE.", call. = FALSE)
    }

    directory <- tdhia_geo_data_dir(geo_accession)

    if (!dir.exists(directory)) {
        message("Directory does not exist: ", directory)
        return(invisible(TRUE))
    }

    if (is.null(geo_accession) && confirm) {
        answer <- readline(
            paste0(
                "Delete all GEO data downloaded by tdhia from\n",
                directory,
                "\n? [y/N]: "
            )
        )

        if (!tolower(trimws(answer)) %in% c("y", "yes")) {
            message("Deletion cancelled.")
            return(invisible(FALSE))
        }
    }

    unlink(
        directory,
        recursive = TRUE,
        force = TRUE
    )

    removed <- !dir.exists(directory)

    if (!removed) {
        warning(
            "Could not completely remove directory: ",
            directory,
            call. = FALSE
        )
    }

    invisible(removed)
}
