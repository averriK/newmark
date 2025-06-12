# nolint start
#' Generate a random sample from a partial-quantile approach, optionally returning fitted parameters.
#'
#' @param n integer, number of samples to generate
#' @param meanValue numeric scalar, the known mean of the distribution in the chosen scale (log-scale, if variable is log-transformed).
#' @param p,q numeric vectors of equal length: the probabilities p and corresponding quantiles q in that same scale.
#' @param parameters logical; if FALSE (default), return only a vector of samples. If TRUE, return a list with (mu, sd, sample).
#' @param nuStart numeric, starting df for the Student-t approach (if used in aggregator).
#' @param deltaD numeric, acceptance band for aggregator methods.
#' @return If \code{parameters=FALSE}, a numeric vector of length \code{n} (the random sample).
#'   If \code{parameters=TRUE}, a list with elements:
#'   \item{muX}{the input meanValue}
#'   \item{sdValue}{the fitted standard deviation from partial quantiles (\code{sdQ()})}
#'   \item{X}{the random sample of length \code{n}}
#' @export rnormQ
rnormQ <- function(
    n = 1,
    meanValue,
    p, q,
    parameters = FALSE,
    nuStart = 10,
    deltaD = 0.01) {
    checkInputs(meanValue, p, q)

    OUT <- fitAllMethodsQ(
        meanValue = meanValue,
        p = p,
        q = q,
        nuStart = nuStart
    )
    agg <- aggregateSigmaQ(OUT, deltaD = deltaD)
    sdValue <- agg$sdValue

    ## --------------------------------------------------------------------
    ## 1)  If the aggregator’s “winner” is method 6 (piecewise CDF), sample
    ##     directly from that piecewise law
    ## --------------------------------------------------------------------
    if (agg$bestMethod == "6") {
        X <- samplePiecewiseY(n, p, q, meanValue, sdValue)
        if (!parameters) {
            return(X)
        }
        return(list(X = X, muX = meanValue, sdValue = sdValue))
    }

    ## --------------------------------------------------------------------
    ## 2)  All other cases proceed exactly as before.
    ## --------------------------------------------------------------------
    if (is.na(sdValue)) {
        ## fall back to pure piecewise if even σ is missing
        sdValue <- OUT$sigma[OUT$method == "6"]
        X <- samplePiecewiseY(n, p, q, meanValue, sdValue)
    } else {
        epsilon <- simulateEps(n, agg$bestMethod, OUT)
        X <- meanValue + sdValue * epsilon
    }

    if (!parameters) {
        return(X)
    }
    list(X = X, muX = meanValue, sdValue = sdValue)
}
