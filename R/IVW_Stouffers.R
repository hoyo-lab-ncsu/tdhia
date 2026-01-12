#' Purpose: To evaluate coordinated methylation effects across CpG within each
#' imprinting control region (ICR) using these steps:
#'  0. Use analyze_associaiton to fit glm model for each CpG. The output
#'  has the Estimate, StdError, Statistic (Walds Z), P_VAL and ADJ_P_VAL
#'  1. Identify "ICRs of interest" (at least one CpG with FDR (ADJ_P_VAL) < 0.1
#'  within the ICR)
#'  2. IVW: inverse-variance-weighted. Typically a meta-analysis method, here
#'  we are treating each CpG as an individual "study"
#'  3. Stouffer's method:
#'  4. Calculate the effect trait (integrates IVW + Stouffer's into one test
#'  statistic)
#'
#' @param aa_output the output dataframe from analyze_association which has
#' fit the lm model for each CpG
#'
#' @returns test
#' @export
#' @author Kate Everly


library(dplyr)
library(tidyverse)

ivw_stouffers <- function(aa_output, verbose=FALSE){
  verbosecat <- function(x) if (verbose) cat(x)
  icr_mapping = tdhia::mapping_cpg_icr_ids # load in map
  cpg_mapping <- icr_mapping %>% # will select just the CpG and ICR ids from the icr_mapping file
    dplyr::select(CpG_id, ICR_id)

  # Step 1: Select ICRs of interest by any ICR that has at least 1 CpG with ADJ_P_VAL<0.10
  cpgs_of_interest <- aa_output %>%
    dplyr::filter(ADJ_P_VAL <0.10) %>%
    dplyr::pull(Variable)
  icr_of_interest <- cpg_mapping %>%
    dplyr::filter(CpG_id %in% cpgs_of_interest) %>%
    dplyr::pull(ICR_id) %>%
    unique()

  # Step 2 and 3: IVW + Stouffer's
  # I originally had a for loop but this apparently uses less memory
  results <- purrr::map_dfr(icr_of_interest, function(icr) {
    # Get the CpG IDs of the CpGs for an ICR
    subset_cpg_ids <- cpg_mapping %>%
      dplyr::filter(ICR_id == icr) %>%
      dplyr::pull(CpG_id)
    #filter aa_output for just subset_cpg_ids
    subset_aa_output <- aa_output %>%
      dplyr::filter(Variable %in% subset_cpg_ids)

    if (nrow(subset_aa_output) == 0) {
      return(tibble::tibble(
        ICR_id = icr,
        Beta_IVW = NA_real_,
        z_icr = NA_real_
      ))
    }
    #summarize the test statistics
    Beta_IVW = sum(subset_aa_output$Estimate / subset_aa_output$StdError^2) / sum(1 / subset_aa_output$StdError^2)
    Z_ICR = sum((1 / subset_aa_output$StdError^2) * subset_aa_output$Statistic) / sqrt(sum((1 / subset_aa_output$StdError^2)^2))
    Effect_trait = -log10(2 * pnorm(-abs(Z_ICR))) * sign(Beta_IVW)

    # Bind results in dataframe
    tibble::tibble(
      ICR_id = icr,
      Beta_IVW = Beta_IVW,
      Z_ICR = Z_ICR,
      Effect_trait = Effect_trait ) })
  return(results)
}

