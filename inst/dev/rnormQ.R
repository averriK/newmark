#' Sample from distribution reconstructed from partial quantiles
#'
#' @param n Number of samples
#' @param meanValue Mean of the distribution
#' @param p Vector of probabilities (0 < p < 1)
#' @param q Vector of quantiles corresponding to p
#' @param parameters Logical, return parameters or just the sample
#' @param nuStart Starting degrees of freedom for t-distribution
#' @param deltaD Acceptance band for sd aggregator
#'
#' @return Either a numeric vector of samples or a list with mean, sd, and sample
#' @export
rnormQ <- function(n = 1,
                   meanValue,
                   p, q,
                   parameters = FALSE,
                   nuStart = 10,
                   deltaD = 0.01) {
  checkInputs(meanValue, p, q)
  OUT <- fitAllMethodsQ(meanValue, p, q, nuStart)
  agg <- aggregateSigmaQ(OUT, deltaD)
  sdValue <- agg$sdValue

  X <- if (agg$bestMethod == "6" || is.na(sdValue)) {
    samplePiecewiseY(n, p, q, meanValue, sdValue)
  } else {
    epsilon <- simulateEps(n, agg$bestMethod, OUT)
    meanValue + sdValue * epsilon
  }

  if (!parameters) return(X)
  list(X = X, muX = meanValue, sdValue = sdValue)
}