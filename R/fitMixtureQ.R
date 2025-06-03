# nolint start
#' Mixture-like aggregator from multiple scenario means plus merged quantiles
#'
#' This helper combines a vector of means (from multiple scenarios, if any)
#' into a single mean via \code{mean(meanValue)}, then fits a distribution
#' to a single set of partial quantiles \code{(p, q)}. The resulting
#' standard deviation is chosen by your existing aggregator logic
#' (see \code{\link{fitAllMethodsQ}} and \code{\link{aggregateSigmaQ}}).
#'
#' @param meanValue numeric vector of means (one per scenario).  Will be averaged.
#' @param p numeric vector of probabilities (strictly between 0 and 1).
#' @param q numeric vector of quantiles corresponding to \code{p}.
#' @param nuStart numeric scalar, starting df for Student-t fit (method 4). Default 10.
#' @param deltaD numeric scalar, acceptance band for aggregator. Default 0.01.
#' @return A single-row \code{data.frame} in the style of \code{\link{templateDF}}:
#'   \item{method}{"mix"}
#'   \item{sigma}{the robust \code{sd} of the final distribution fit}
#'   \item{d}{\code{NA}}
#'   \item{sse}{\code{NA}}
#'   \item{muHat}{the final mean (average of \code{meanValue})}
#'   \item{other columns}{\code{NA}}
#'
#' @keywords internal
#' @noRd
fitMixtureQ <- function(
    meanValue,
    p,
    q,
    nuStart = 10,
    deltaD  = 0.01
) {
  # 1) Validate user input
  if (!is.numeric(meanValue) || !all(is.finite(meanValue))) {
    stop("`meanValue` must be numeric and finite.")
  }
  if (length(p) != length(q)) {
    stop("`p` and `q` must have the same length.")
  }
  if (any(p <= 0 | p >= 1)) {
    stop("All probabilities in `p` must lie strictly between 0 and 1.")
  }
  
  # Ensure p, q are sorted in ascending order of p
  if (is.unsorted(p)) {
    o <- order(p)
    p <- p[o]
    q <- q[o]
  }
  
  # 2) Aggregate the multiple means (simple average)
  muComb <- mean(meanValue)
  
  # 3) Re-use your standard pipeline:
  #    fitAllMethodsQ => aggregateSigmaQ => fallback if needed
  OUT <- fitAllMethodsQ(
    meanValue = muComb,
    p         = p,
    q         = q,
    nuStart   = nuStart
  )
  agg <- aggregateSigmaQ(OUT, deltaD)
  sdValue <- agg$sdValue
  # fallback check typically done inside aggregateSigmaQ
  
  # 4) Build the final single-row data.frame
  #    method="mix", store the combined mean in muHat, set d/sse=NA
  out <- templateDF(
    method = "mix",
    sigma  = sdValue,
    d      = NA_real_,
    sse    = NA_real_,
    muHat  = muComb
  )
  
  return(out)
}
# nolint end
