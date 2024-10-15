#' Alzheimer’s Disease Brain Methylation Data
#'
#' This dataset contains the probe-level beta values for 16 Alzheimer’s Disease brain samples.
#'
#' @name ad_brain_n16_probe_beta_matrix
#' @docType data
#' @keywords datasets
#' @usage NULL
#' @format A matrix with rows as probes and columns as samples. The columns represent the samples, and each column contains raw beta values for that sample.
#' \describe{
#'   \item{Column names}{Sample IDs corresponding to the IDAT files.}
#'   \item{Metadata}{
#'     The following table describes the disease status and race of each sample:
#'     \tabular{lll}{
#'       \strong{Sample Name} \tab \strong{Disease} \tab \strong{Race} \cr
#'       207344530004_R01C01 \tab Control \tab black \cr
#'       207344530004_R03C01 \tab Control \tab black \cr
#'       207344530004_R05C01 \tab Alzheimer_Disease \tab black \cr
#'       207344530004_R07C01 \tab Control \tab black \cr
#'       207344530004_R09C01 \tab Alzheimer_Disease \tab black \cr
#'       207344530004_R11C01 \tab Control \tab black \cr
#'       207344530004_R02C01 \tab Alzheimer_Disease \tab black \cr
#'       207344530004_R04C01 \tab Alzheimer_Disease \tab black \cr
#'       207344530004_R06C01 \tab Control \tab white \cr
#'       207344530004_R08C01 \tab Alzheimer_Disease \tab white \cr
#'       207344530004_R10C01 \tab Control \tab white \cr
#'       207344530004_R12C01 \tab Control \tab white \cr
#'       207344530004_R01C02 \tab Alzheimer_Disease \tab white \cr
#'       207344530004_R03C02 \tab Control \tab white \cr
#'       207344530004_R05C02 \tab Alzheimer_Disease \tab white \cr
#'       207344530004_R07C02 \tab Alzheimer_Disease \tab white
#'     }
#'   }
#' }
#' @source 16 Alzheimer's brain samples with matching WGBS that used in array design paper
"ad_brain_n16_probe_beta_matrix"

