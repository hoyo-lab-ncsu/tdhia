

#' @title False Discovery Rate Correction with Custom Number of Comparisons 
#'
#' @description
#' Applies the Benjamini-Hochberg (BH) procedure for False Discovery Rate (FDR)
#' correction, using a user-specified value for the total number of comparisons (n),
#' which may be less than the actual number of p-values provided.
#'
#' @param p_values A numeric vector of raw p-values.
#' @param n A positive integer specifying the total number of
#' comparisons (n) to use in the BH formula.
#'
#' @return A numeric vector of adjusted p-values (q-values), ordered to
#' correspond to the input p_values.
#' @export
apply_fdr_correction <- function(pvals, n) {
  # pvals: vector of raw p-values
  # n: number of independent tests to assume (≤ length(pvals))
  # method: any method accepted by p.adjust (e.g., "BH", "bonferroni", "holm", etc.)
  
  stopifnot(length(pvals) >= n,
            all(pvals >= 0 & pvals <= 1))
  
  # Adjust p-values as if only n tests were performed
  
  # Benjamini-Hochberg adjustment with custom number of comparisons
  ord <- order(pvals)
  ranks <- rank(pvals, ties.method = "first")
  adj_p <- pvals * n / ranks
  adj_p <- pmin(adj_p, 1)
  adj_p[ord] <- cummin(adj_p[ord[length(ord):1]])[length(ord):1]  # monotonic correction
  
  
  return(adj_p)
}




#' custom_p.adjust
#'
#' @description Given a set of p-values, returns p-values adjusted using one of
#' several methods.
#'
#' @param p numeric vector of p-values (possibly with NAs). Any other R object is coerced by as.numeric.
#' @param method correction method, a character string. Can be abbreviated.
#' Options include "bonferroni", "holm", "hochberg", "hommel", "BH", "BY", or "none".
#' @param n of comparisons.
#'
#' @importFrom stats p.adjust
#' @export
custom_p.adjust <- function (p, method = stats::p.adjust.methods, n = length(p))
{
  method <- match.arg(method)
  if (method == "fdr")
    method <- "BH"
  nm <- names(p)
  p <- as.numeric(p)
  p0 <- stats::setNames(p, nm)
  if (all(nna <- !is.na(p)))
    nna <- TRUE
  else p <- p[nna]
  lp <- length(p)
  # Commented this out to allow reduced n with p.value adjustment
  # See: https://stackoverflow.com/questions/30108510/p-adjust-with-n-than-number-of-tests
  # stopifnot(n >= lp)
  if (n <= 1)
    return(p0)

  if (n == 2 && method == "hommel")
    method <- "hochberg"

  p0[nna] <- switch(method, bonferroni = pmin(1, n * p), holm = {
    i <- seq_len(lp)
    o <- order(p)
    ro <- order(o)
    pmin(1, cummax((n + 1L - i) * p[o]))[ro]
  }, hommel = {
    if (n > lp) p <- c(p, rep.int(1, n - lp))
    i <- seq_len(n)
    o <- order(p)
    p <- p[o]
    ro <- order(o)
    q <- pa <- rep.int(min(n * p/i), n)
    for (j in (n - 1L):2L) {
      ij <- seq_len(n - j + 1L)
      i2 <- (n - j + 2L):n
      q1 <- min(j * p[i2]/(2L:j))
      q[ij] <- pmin(j * p[ij], q1)
      q[i2] <- q[n - j + 1L]
      pa <- pmax(pa, q)
    }
    pmax(pa, p)[if (lp < n) ro[1L:lp] else ro]
  }, hochberg = {
    i <- lp:1L
    o <- order(p, decreasing = TRUE)
    ro <- order(o)
    pmin(1, cummin((n + 1L - i) * p[o]))[ro]
  }, BH = {
    i <- lp:1L
    o <- order(p, decreasing = TRUE)
    ro <- order(o)
    pmin(1, cummin(n/i * p[o]))[ro]
  }, BY = {
    i <- lp:1L
    o <- order(p, decreasing = TRUE)
    ro <- order(o)
    q <- sum(1/(1L:n))
    pmin(1, cummin(q * n/i * p[o]))[ro]
  }, none = p)
  # Added this line to prevent adjusted p-values that are smaller than unadjusted.
  p0 <- ifelse(p0 < p, p, p0)
  p0
}
