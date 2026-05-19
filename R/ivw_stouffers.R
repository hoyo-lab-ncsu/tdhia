#' ivw_stouffers
#'
#' Purpose: To evaluate coordinated methylation effects across CpG within each
#' imprinting control region (ICR) using these steps:
#'  0. Use analyze_associaiton to fit glm model for each CpG. The output
#'  has the Variable (CpG), Estimate, StdError, Statistic (Walds Z), P_VAL
#'  and ADJ_P_VAL
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
#' @param adj_p_level numeric. Max ADJ_P_VAL for a CpG to be included in CpGs/ICRs of
#' interest. Default = 0.10
#' @param min_cpg minimum number of cpg sites for an icr to be included in
#' results (default is 3).
#' @param verbose boolean, when TRUE, prints progress to command line (default = TRUE).
#' @param top_n_icrs number of top ICRs you would like the final results to be
#' filtered to include. Default=NULL (all ICRs reported)
#'
#' @returns test
#' @export
#' @author Kate Everly


ivw_stouffers <- function(aa_output, adj_p_level=0.10, min_cpg=3, verbose=TRUE, top_n_icrs=NULL){
  verbosecat <- function(x) if (verbose) cat(x)
  icr_mapping = tdhia::mapping_cpg_icr_ids # load in map
  cpg_mapping <- icr_mapping %>% # will select just the CpG and ICR ids from the icr_mapping file
    dplyr::select(CpG_id, ICR_id) %>%
    dplyr::distinct()

  # Extra safety steps from Dereje's code:
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # ----Safety checks & helpers ----
  required_cols <- c("Variable","Estimate","StdError","Statistic","P_VAL","ADJ_P_VAL")
  stopifnot(all(required_cols %in% colnames(aa_output)))

  # If CpG-level Z is missing/non-finite, recompute = Estimate / StdError
  if (!"Statistic" %in% names(aa_output) || any(!is.finite(aa_output$Statistic))) {
    message("Recomputing CpG-level Z as Estimate/StdError where needed.")
    aa_output$Statistic <- with(aa_output,
                                ifelse(is.finite(Estimate) & is.finite(StdError) & StdError > 0,
                                Estimate / StdError, NA_real_)
      )
    }
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


  # Step 1: Select ICRs of interest by any ICR that has at least 1 CpG with
  # ADJ_P_VAL <0.10 (or chosen adj_p_val)
  cpgs_of_interest <- aa_output %>%
    dplyr::filter(ADJ_P_VAL < adj_p_level) %>%
    dplyr::pull(Variable)

  icr_of_interest <- cpg_mapping %>%
    dplyr::filter(CpG_id %in% cpgs_of_interest) %>%
    dplyr::pull(ICR_id) %>%
    unique()

  # Step 2 and 3: IVW + Stouffer's
  # I originally had a for loop but purrr apparently uses less memory
  results <- purrr::map_dfr(icr_of_interest, function(icr) {
    # Get the CpG IDs of the CpGs for an ICR
    subset_cpg_ids <- cpg_mapping %>%
      dplyr::filter(ICR_id == icr, CpG_id %in% aa_output$Variable) %>%
      dplyr::pull(CpG_id) %>%
      unique()

    #filter aa_output for just subset_cpg_ids
    subset_aa_output <- aa_output %>%
      dplyr::filter(Variable %in% subset_cpg_ids)

    num_cpgs_in_icr <- nrow(subset_aa_output)

    #the subset aa_output should never have 0 rows bc. you already filtered for sig CpGs....
    # just leaving in anyways
    if (nrow(subset_aa_output) == 0) {
      return(tibble::tibble(
        ICR_id = icr,
        n_CpGs = 0,
        Beta_IVW = NA_real_,
        Z_ICR = NA_real_,
        ICR_P_VAL = NA_real_,
        Effect_trait = NA_real_
      ))
    }

    # Summarize the test statistics and return NA if result is infinite or zero
    weights <- 1 / subset_aa_output$StdError^2
    Beta_IVW <- sum(weights * subset_aa_output$Estimate) / sum(weights)
    Z_ICR <- sum(weights * subset_aa_output$Statistic) / sqrt(sum(weights^2))

    p_ICR <- 2 * pnorm(-abs(Z_ICR))
    Effect_trait <- -log10(p_ICR) * sign(Beta_IVW)

    # Validity check- if there's values of infinity or 0 (or negative for p-value),
    # that ICR will just get NAs
    if (!is.finite(Beta_IVW) ||
        !is.finite(Z_ICR)    ||
        !is.finite(Effect_trait) ||
        !is.finite(p_ICR) ||
        p_ICR <= 0) {
      return(tibble::tibble(
        ICR_id = icr,
        n_CpGs = num_cpgs_in_icr,
        Beta_IVW = NA_real_,
        Z_ICR = NA_real_,
        ICR_P_VAL = NA_real_,
        Effect_trait = NA_real_
      ))
    }

    # Bind results in dataframe
    tibble::tibble(
      ICR_id = icr,
      n_CpGs = num_cpgs_in_icr,
      Beta_IVW = Beta_IVW,
      Z_ICR = Z_ICR,
      ICR_P_VAL = p_ICR,
      Effect_trait = Effect_trait
    )
  })
  # Filter for only ICRs with at least the min number of CpGs
  verbosecat(sprintf("> Filtering out ICRs with < %d cpg sites...\n", min_cpg))
  results <- results %>%
    dplyr::filter(n_CpGs >= min_cpg)

  if (is.null(top_n_icrs)) {
    return(results)
  }

  if (!is.null(top_n_icrs)) {
    results_ranked <- results %>%
      filter(is.finite(Effect_trait ), is.finite(ICR_P_VAL)) %>%
      mutate(Abs_Effect_trait = abs(Effect_trait)) %>%
      arrange(desc(Abs_Effect_trait), ICR_P_VAL, ICR_id) %>%
      distinct(ICR_id, .keep_all = TRUE)
    # Take exactly topN (no overshoot on ties)
    top_icr_results <- results_ranked %>%
      slice_head(n = top_n_icrs) #%>% pull(ICR_id)
    return(top_icr_results)
  }
}
