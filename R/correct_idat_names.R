
#' correct_idat_names
#'
#' @description correct capitalization to match the convention expected by the
#' sesame package. Example correct filenames:
#'   4207113116_1544_Grn.idat.
#'   4207113116_1544_Red.idat.
#' @param idat_dir_path either a full filesystem path to the directory of idats
#' to be rename, or a vector paths to individual IDAT files
#' @param rename boolean, when true IDAT files are renamed, otherwise their
#' fileneames are checked and an error or warning is raised
#' @param failcheck_error boolean, when true an error is raised when filenames are
#' checked, otherwise a warning is raised.
#'
#' @return results with fitted coefficients from the glm, sorted by p-value
correct_idat_names <- function (idat_dir_path, rename = TRUE, failcheck_error = TRUE) {

  if (length(idat_dir_path)==1) {
    # Get all idat files in target directory
    idat_fullnames <- dir(path = idat_dir_path, pattern = "/*.idat",
                          full.names = TRUE, ignore.case = TRUE)
  } else {
    idat_fullnames <- idat_dir_path
  }


  if (!rename) {
    pass_bv <- grepl("_(Grn)|(Red).idat", x = idat_fullnames, ignore.case = FALSE )
    if (all(pass_bv)) return(TRUE)

    cat("These IDAT files do not match proper naming convention:\n")
    cat(idat_fullnames[!pass_bv], sep = "\n")

    if (failcheck_error){
      stop("correct_idat_names: IDAT files not properly named, rename with correct_idat_names()")
    } else {
      warning("correct_idat_names: IDAT files not properly named, rename with correct_idat_names()")
      return(FALSE)
    }

  } else {
    cat("Renaming IDAT Files...\n")
    # For IDAT extention lower case
    new_idat_fullnames = stringr::str_replace(idat_fullnames, '.IDAT$', '.idat')
    # Force green camel case
    new_idat_fullnames = stringr::str_replace(new_idat_fullnames, 'grn.idat', 'Grn.idat')
    # Force red camel case
    new_idat_fullnames = stringr::str_replace(new_idat_fullnames, 'red.idat', 'Red.idat')

    # organize into dataframe
    df_log = data.frame(old_name = idat_fullnames, new_name = new_idat_fullnames)

    # Write log of proposed renames
    utils::write.csv(x=df_log, file = paste0(idat_dir_path, '/idat_rename_log.csv'))

    # Rename files
    file.rename(from = df_log$old_name, to = df_log$new_name)
  }

}



