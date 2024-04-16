#' simulate_idat_data
#'
#' Produces a series of simulated sigsets (object representing imported IDAT file
#' data in SeSame package). Uses inter
#' 
#' simulate_idat_data().
#' 
#' Data is stored as a simple ecfd, later versions should switch to a smoothed 
#' pdf of the fluorescent values.
#'
#' @param n_idats number of sigsets to simulate
#' @param sigset_model model of idat file to be used as reference to simulations.
#' @param rand.seed random seed to reproduce results
#' @param data_colnames column names for fluorescent data
#' @return sigset list: list of sigset objects that can be used as input to 
#'    openSesame().
#' 
simulate_idat_data <- function(n_idats=1, sigset_model = NULL, rand.seed = 0,
                               data_colnames = c("MG", "MR", "UG", "UR")) {

  # If no sigset_model is specified, loads internal one
  if (is.null(sigset_model)) { sigset_model = tdhia::sigset_model}

  # Set random seed for reproducibility
  set.seed(rand.seed)
  
  # Use sigset model to sample fluorescent data, data column by data column
  # Then set missing data in original ref file to NA in simulated files.
  sigset_list = list()
  for (k in 1:n_idats) {
    # Initially copy sigset, then replace data
    sigset_list[[k]] <- sigset_model$sigset_template
    
    for (c in seq_along(data_colnames)) {
       temp <- 
        sample(sigset_model$sigset_values[[data_colnames[c]]], 
               size = nrow(sigset_model$sigset_template), replace = TRUE)
       temp[sigset_model$sigset_mask[[data_colnames[c]]]] <- NA
    sigset_list[[k]][[data_colnames[c]]] <- temp
    }
    
  }

  
  return(sigset_list)
}
