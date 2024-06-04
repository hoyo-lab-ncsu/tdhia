



#' copy_idat_files
#' 
#' @description
#'  Given a input folder path, and a list of idat basenames, copys the red and
#'  green IDAT files to the export file direcotry. Useful for subsetting the 
#'  IDAT files.
#'  
#'  @param idat_basenames list of basenames for the IDAT files to be copied
#'  @param idat_dir_paths full system path to IDAT files
#'  @param export_file_dir output path to copy IDAT files.
#' 
#' 
copy_idat_files <- function(idat_basenames, idat_dir_paths, export_file_dir) {


  # Get all IDAT files associated with idat basenames
  idate_filenames_list <- base::sapply(X = idat_basenames,FUN = function(x)
    base::list.files(path = idat_dir_paths, pattern = paste0(x, '.*.idat$'),
                     ignore.case = TRUE))
  idate_filenames_list <- base::as.vector((idate_filenames_list))
  
  # Copy files
  res <- file.copy(from = paste0(idat_dir_paths, "/",idate_filenames_list),
            to = export_file_dir, overwrite = TRUE)
  
 return(res)
}