#' Plot density of beta values for a subset of samples
#'
#' @param subdata A data frame containing beta values, with probes as rows and samples as columns.
#' @param remove_na A logical value indicating whether to remove NA values.
#' @param show_legend A logical value indicating whether to show the legend.
#' @param legend_position A string indicating the position of the legend. Options are "none", "top", "bottom", "left", and "right".
#' @return A ggplot object representing the density plot.
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom stats na.omit
#' @export
#' @examples
#' # Load required namespaces
#' requireNamespace("ggplot2", quietly = TRUE)
#' requireNamespace("reshape2", quietly = TRUE)
#' # Assuming 'subdata' is your data frame
#' # Create example data frame with row names as probeset_id (replace this with your actual data)
#' subdata <- data.frame(
#'   Sample1 = runif(100),
#'   Sample2 = runif(100),
#'   Sample3 = runif(100)
#' )
#' rownames(subdata) <- paste0("probe", 1:100)  # Replace "probe" with your actual prefix
#'
#' # Let's say we want to plot the density and remove NA values,
#' # and show the legend on the bottom
#' density_plot <- plot_beta_density(subdata,
#'                                  remove_na = TRUE,
#'                                  show_legend = TRUE,
#'                                  legend_position = "right")
#' print(density_plot)
plot_beta_density <- function(subdata, remove_na = TRUE, show_legend = TRUE, legend_position = "right") {

    # Melt the data frame
  rawbeta <- melt(as.matrix(subdata))

  # Rename the variable column
  colnames(rawbeta)[2] <- "SampleName"

  # Remove NA values if specified
  if (remove_na) {
    rawbeta <- na.omit(rawbeta)
  }

  # Determine legend options
  legend_options <- if (show_legend) {
    list(legend.position = legend_position)
  } else {
    list(legend.position = "none")
  }

  # Plot the density of beta values for each sample with theme_minimal() and legend options
  p <- ggplot2::ggplot(rawbeta, aes_string(x = "value", color = "SampleName")) +
    geom_density() +
    labs(x = "Beta Value", y = "Density") +
    theme_minimal() +
    theme(legend.position = legend_options$legend.position)

  return(p)
}
