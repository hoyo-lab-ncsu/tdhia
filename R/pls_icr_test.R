

#' PLS-based CpG-set association test with covariates and resampling
#'
#' Performs a Partial Least Squares (PLS) CpG-set test for a single ICR/region,
#' adjusting for sample-level covariates and providing valid inference via
#' resampling. Continuous outcomes use Freedman–Lane residual permutation;
#' binary outcomes use a parametric bootstrap under the null logistic model.
#'
#' @param y Outcome vector of length n.
#'   For \code{outcome = "continuous"}, numeric.
#'   For \code{outcome = "binary"}, must be 0/1, logical, or a 2-level factor.
#' @param Z Numeric matrix of CpG methylation values for a single ICR
#'   with dimensions n (samples) x p (CpGs).
#' @param X Optional data.frame or matrix of sample-level covariates
#'   with n rows (e.g., age, sex, BMI, batch).
#'   If NULL, an intercept-only model is used.
#' @param k Integer number of PLS components to retain (latent dimensions).
#'   Must be pre-specified. Typical values are 1 or 2.
#' @param B Integer number of resampling iterations.
#'   For continuous outcomes, number of permutations.
#'   For binary outcomes, number of bootstrap replicates.
#' @param outcome Character string specifying outcome type.
#'   Either \code{"continuous"} or \code{"binary"}.
#' @param scale_Z Logical; whether to scale CpGs inside the PLS step.
#' @param center_Z Logical; whether to center CpGs inside the PLS step.
#' @param seed Integer random seed for reproducibility.
#' @param return_fits Logical; if TRUE, returns fitted model objects
#'   (null model, alternative model, and PLS fit).
#' @return A list containing test statistics, resampling p-values,
#' PLS coefficients, and component diagnostics.
#'
#' @importFrom stats lm glm anova fitted residuals coef deviance pchisq rbinom
#'
#' @export
pls_icr_test <- function(y, Z, X = NULL,
    k = 1, B = 2000, outcome = c("continuous", "binary"),
    scale_Z = TRUE, center_Z = TRUE, seed = 1,
    return_fits = FALSE, db_flag = F) 
  {
  
  # ------------------------------------------------------------
  # 0. Basic argument checks and setup
  # ------------------------------------------------------------
  
  outcome <- match.arg(outcome)
  
  # Ensure CpG matrix is numeric and matrix-like
  Z <- as.matrix(Z)
  n <- nrow(Z)
  
  # Outcome length must match sample count
  if (length(y) != n) {
    stop("Length of y must match nrow(Z).")
  }
  
  # Detect whether covariates are present
  has_cov <- !is.null(X)
  
  # If covariates exist, coerce to data.frame and check dimensions
  if (has_cov) {
    Xdf <- as.data.frame(X)
    if (nrow(Xdf) != n) {
      stop("nrow(X) must match nrow(Z).")
    }
  }
  
  # ------------------------------------------------------------
  # 1. Outcome preprocessing + missingness handling
  # ------------------------------------------------------------
  # We construct a complete-case index that simultaneously removes
  # missing values from:
  #   - outcome y
  #   - CpG methylation matrix Z
  #   - covariates X (if present)
  #
  # This ensures all downstream models are fit on identical samples.
  
  if (outcome == "continuous") {
    # Continuous outcome must be numeric
    y_vec <- as.numeric(y)
    # Complete-case filtering across y, Z, and X
    cc <- stats::complete.cases(y_vec, Z, if (has_cov) Xdf else rep(TRUE, n))
    
    y_vec <- y_vec[cc]
    
  } else {
    
    # Binary outcome: accept 0/1 numeric, logical, or 2-level factor
    y_bin <- y
    
    if (is.factor(y_bin) || is.character(y_bin)) {
      y_bin <- as.integer(as.factor(y_bin)) - 1L
    } else if (is.logical(y_bin)) {
      y_bin <- as.integer(y_bin)
    } else {
      y_bin <- as.numeric(y_bin)
    }
    
    # Enforce valid Bernoulli support
    ok01 <- is.finite(y_bin) & y_bin %in% c(0, 1)
    
    cc <- stats::complete.cases(
      y_bin,
      Z,
      if (has_cov) Xdf else rep(TRUE, n)
    ) & ok01
    
    y_vec <- y_bin[cc]
  }
  
  # Apply the same filtering to Z and X
  Z <- Z[cc, , drop = FALSE]
  if (has_cov) Xdf <- Xdf[cc, , drop = FALSE]
  n <- nrow(Z)
  
  # Guard against degenerate cases
  if (n < k + 5) { stop("Too few samples relative to k.")  }
  
  set.seed(seed)
  
  # ------------------------------------------------------------
  # 2. Helper: fit PLS and extract k latent scores
  # ------------------------------------------------------------
  # This function:
  #   - runs PLS of a pseudo-response r on CpGs Z
  #   - extracts the first k latent dimensions (scores)
  #
  # IMPORTANT: This step *uses outcome information* (via r),
  # which is why it must be recomputed inside resampling.
  
  fit_pls <- function(r) {
    
    pls_fit <- pls::plsr(
      r ~ Z,
      ncomp = k,
      method = "simpls",
      center = center_Z,
      scale = scale_Z,
      validation = "none"
    )
    
    # Extract latent scores (n x k)
    scores <- pls_fit$scores[, seq_len(k), drop = FALSE]
    colnames(scores) <- paste0("PLS", seq_len(k))
    
    list(scores = scores, fit = pls_fit)
  }
  
  # ============================================================
  # 3. Continuous outcome case (Gaussian LM)
  # ============================================================
  if (outcome == "continuous") {
    
    # ----------------------------------------------------------
    # 3a. Null model: y ~ X
    # ----------------------------------------------------------
    # This removes covariate effects before PLS is applied.
    # Freedman–Lane permutation requires residuals from this model.
    
    null_fit <- if (has_cov) {
      stats::lm(y_vec ~ ., data = cbind.data.frame(y = y_vec, Xdf))
    } else {
      stats::lm(y_vec ~ 1)
    }
    
    y_hat <- stats::fitted(null_fit)     # fitted values under null
    e_hat <- stats::residuals(null_fit)  # covariate-adjusted residuals
    
    # ----------------------------------------------------------
    # 3b. Observed PLS on residualized outcome
    # ----------------------------------------------------------
    # This ensures PLS components represent association *beyond*
    # what covariates explain.
    
    pls_obs <- fit_pls(e_hat)
    
    # ----------------------------------------------------------
    # 3c. Alternative model: y ~ X + PLS scores
    # ----------------------------------------------------------
    alt_fit <- stats::lm(
      y_vec ~ .,
      data = cbind.data.frame(y = y_vec, Xdf, pls_obs$scores)
    )
    
    # Partial F-statistic for adding k PLS dimensions
    aobs <- stats::anova(null_fit, alt_fit)
    F_obs <- unname(aobs$F[2])
    
    # Naive p-value (INVALID for inference; diagnostic only)
    p_naive <- unname(aobs$`Pr(>F)`[2])
    
    # ----------------------------------------------------------
    # 3d. Freedman–Lane permutation
    # ----------------------------------------------------------
    # Permute null residuals, reconstruct y*, and re-run:
    #   null -> PLS -> alternative
    # This preserves covariate structure and CpG correlation.
    
    F_perm <- numeric(B)
    
    for (b in seq_len(B)) {
      
      # Permute residuals
      e_perm <- sample(e_hat)
      
      # Reconstruct permuted outcome
      y_perm <- y_hat + e_perm
      
      # Refit null on permuted outcome
      null_b <- if (has_cov) {
        stats::lm(y_perm ~ ., data = cbind.data.frame(y_perm, Xdf))
      } else {
        stats::lm(y_perm ~ 1)
      }
      
      # PLS must be recomputed under permutation
      pls_b <- fit_pls(stats::residuals(null_b))
      
      # Fit alternative model
      alt_b <- stats::lm(
        y_perm ~ .,
        data = cbind.data.frame(y_perm, Xdf, pls_b$scores)
      )
      
      # Record partial F-statistic
      ab <- stats::anova(null_b, alt_b)
      F_perm[b] <- unname(ab$F[2])
    }
    
    # Empirical p-value
    p_resample <- (1 + sum(F_perm >= F_obs)) / (B + 1)
    
    # ----------------------------------------------------------
    # 3e. Directional summary
    # ----------------------------------------------------------
    # Direction is taken from the sign of the first PLS coefficient.
    
    coef_alt <- stats::coef(alt_fit)
    pls_coef <- coef_alt[grep("^PLS", names(coef_alt))]
    direction <- sign(pls_coef[1])
    
    out <- list(
      outcome = "continuous",
      p_resample = p_resample,
      stat_obs = F_obs,
      p_naive = p_naive,
      pls_coef = pls_coef,
      direction_first_pls = direction,
      pls_scores = pls_obs$scores,
      pls_weights = pls_obs$fit$loading.weights[, seq_len(k), drop = FALSE],
      pls_loadings = pls_obs$fit$loadings[, seq_len(k), drop = FALSE]
    )
    
  } else {
    
    # ============================================================
    # 4. Binary outcome case (logistic regression)
    # ============================================================
    # Residual permutation is not valid for Bernoulli outcomes.
    # We instead use a parametric bootstrap under the null.
    
    # ----------------------------------------------------------
    # 4a. Null logistic model: y ~ X
    # ----------------------------------------------------------
    null_fit <- if (has_cov) {
      stats::glm(
        y_vec ~ .,
        data = cbind.data.frame(y = y_vec, Xdf),
        family = stats::binomial()
      )
    } else {
      stats::glm(y_vec ~ 1, family = stats::binomial())
    }
    
    # Fitted probabilities under null
    p_hat <- stats::fitted(null_fit)
    
    # Pseudo-response for PLS (score-like residual)
    r_obs <- y_vec - p_hat
    
    # ----------------------------------------------------------
    # 4b. Observed PLS + alternative model
    # ----------------------------------------------------------
    pls_obs <- fit_pls(r_obs)
    
    alt_fit <- stats::glm(
      y_vec ~ .,
      data = cbind.data.frame(y = y_vec, Xdf, pls_obs$scores),
      family = stats::binomial()
    )
    
    # Likelihood ratio statistic
    LRT_obs <- stats::deviance(null_fit) - stats::deviance(alt_fit)
    
    # Naive chi-square p-value (INVALID; diagnostic only)
    p_naive <- stats::pchisq(LRT_obs, df = k, lower.tail = FALSE)
    
    # ----------------------------------------------------------
    # 4c. Parametric bootstrap under null
    # ----------------------------------------------------------
    LRT_boot <- numeric(B)
    
    for (b in seq_len(B)) {
      
      # Simulate outcome from null model
      y_b <- stats::rbinom(n, 1, p_hat)
      
      # Refit null to simulated data
      null_b <- if (has_cov) {
        stats::glm(
          y_b ~ .,
          data = cbind.data.frame(y_b, Xdf),
          family = stats::binomial()
        )
      } else {
        stats::glm(y_b ~ 1, family = stats::binomial())
      }
      
      # Recompute PLS under bootstrap
      pls_b <- fit_pls(y_b - stats::fitted(null_b))
      
      # Fit alternative model
      alt_b <- stats::glm(
        y_b ~ .,
        data = cbind.data.frame(y_b, Xdf, pls_b$scores),
        family = stats::binomial()
      )
      
      # Record LRT
      LRT_boot[b] <- stats::deviance(null_b) - stats::deviance(alt_b)
    }
    
    # Empirical bootstrap p-value
    p_resample <- (1 + sum(LRT_boot >= LRT_obs)) / (B + 1)
    
    coef_alt <- stats::coef(alt_fit)
    pls_coef <- coef_alt[grep("^PLS", names(coef_alt))]
    direction <- sign(pls_coef[1])
    
    out <- list(
      outcome = "binary",
      p_resample = p_resample,
      stat_obs = LRT_obs,
      p_naive = p_naive,
      pls_coef = pls_coef,
      direction_first_pls = direction,
      pls_scores = pls_obs$scores,
      pls_weights = pls_obs$fit$loading.weights[, seq_len(k), drop = FALSE],
      pls_loadings = pls_obs$fit$loadings[, seq_len(k), drop = FALSE]
    )
  }
  
  # Optionally return fitted model objects for diagnostics
  if (return_fits) {
    out$null_fit <- null_fit
    out$alt_fit  <- alt_fit
    out$pls_fit  <- pls_obs$fit
  }
  
  out
}
