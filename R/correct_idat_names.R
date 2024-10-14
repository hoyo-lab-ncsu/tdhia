
#' correct_idat_names
#'
#' @description correct capitalization to match the convention expected by the
#' sesame package. Example filename: 4207113116_A_Grn.idat.
#'
#' @param idat_dir_path full filesystem path to location of idats to be renamed.
#'
#' @return results with fitted coefficients from the glm, sorted by p-value
correct_idat_names <- function (idat_dir_path) {

  # Get all idat files in target directory
  idat_fullnames <- dir(path = idat_dir_path, pattern = "/*.idat",
                            full.names = TRUE, ignore.case = TRUE)

  # For IDAT extention lower case
  new_idat_fullnames = stringr::str_replace(idat_fullnames, '.IDAT$', '.idat')
  # Force green camel case
  new_idat_fullnames = stringr::str_replace(new_idat_fullnames, 'grn.idat', 'Grn.idat')
  # Force red camel case
  new_idat_fullnames = stringr::str_replace(new_idat_fullnames, 'red.idat', 'Red.idat')

  # organize into dataframe
  df_log = data.frame(old_name = idat_fullnames, new_name = new_idat_fullnames)

  # Write log of proposed renames
  write.csv(x=df_log, file = paste0(idat_dir_path, '/idat_rename_log.csv'))

  # Rename files
  file.rename(from = df_log$old_name, to = df_log$new_name)

}



