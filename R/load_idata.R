
#' load_idata
#'
#' \code{load_idata} returns the sum of all the values present in its arguments.
#'
#' This is a generic function: methods can be defined for it directly
#' or via the \code{\link{Summary}} group generic. For this to work properly,
#' the arguments \code{...} should be unnamed, and dispatch is on the
#' first argument.
#' @param idata load the idata from red and green channel
#'
#' @return a matrix samples to load
#' @export
#'
#' @examples
#' load-data(input_data)
load_idata <- function(idata) {
    # write this later
    warning("You ran foo and I haven't written it yet")

}


library(ggplot2)
show_mtcars <- function() {
  qplot(mpg, wt, data = mtcars)
}

# You need the suggested package for this function
my_fun1 <- function(a, b) {
  if (!requireNamespace("pkg", quietly = TRUE)) {
    stop("Pkg needed for this function to work. Please install it.",
         call. = FALSE)
  }
}
# There's a fallback method if the package isn't available
my_fun2 <- function(a, b) {
  if (requireNamespace("pkg", quietly = TRUE)) {
    pkg::f()
  } else {
    g()
  }
}
