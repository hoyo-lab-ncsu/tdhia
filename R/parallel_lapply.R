#' parallel_lapply
#'
#' @description a wrapper to emulate mclapply behavior on Windows systems by using
#' parLapply.
#'
#' @param X a vector (atomic or list) or an expressions vector. Other objects
#' (including classed objects) will be coerced by as.list.
#' @param FUN the function to be applied to each element of X.
#' @param mc.cores he number of cores to use.
#'
#' @return a list of the same length as X and named by X
#'
parallel_lapply <- function(X, FUN, mc.cores,...) {


  if (.Platform$OS.type == "windows" && mc.cores > 1) {
    cat("* Simulating mcapply on windows with parLapply...\n")
    cl <- parallel::getDefaultCluster()
    if (is.null(cl)){
      cl <- parallel::makeCluster( min(mc.cores, parallel::detectCores()) )
    }

    # Obtain list of user defined functions from this repository and export to
    # each cluster
    results1 <- purrr::map_dfr(ls(envir=.GlobalEnv), ~data.frame(name=.,
                                                                 is_func = is.function(get(.))))
    results2 <- filter(results1,is_func)
    parallel::clusterEvalQ(cl, library("dplyr"))
    parallel::clusterExport(cl = cl, varlist = results2$name)

    output <- parallel::parLapply(cl=cl,X=X, fun=FUN,...)

  } else  {
    if (mc.cores == 1) {
      cat("* Executing singlecore mcapply on windows...\n")
    }
    output <- parallel::mclapply(X = X, FUN = FUN, mc.cores = mc.cores,...)
  }

  return(output)
}
