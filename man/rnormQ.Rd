% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnormQ.R
\name{rnormQ}
\alias{rnormQ}
\title{Generate a random sample from a partial-quantile approach, optionally returning fitted parameters.}
\usage{
rnormQ(n = 1, meanValue, p, q, parameters = FALSE, nuStart = 10, deltaD = 0.01)
}
\arguments{
\item{n}{integer, number of samples to generate}

\item{meanValue}{numeric scalar, the known mean of the distribution in the chosen scale (log-scale, if variable is log-transformed).}

\item{p, q}{numeric vectors of equal length: the probabilities p and corresponding quantiles q in that same scale.}

\item{parameters}{logical; if FALSE (default), return only a vector of samples. If TRUE, return a list with (mu, sd, sample).}

\item{nuStart}{numeric, starting df for the Student-t approach (if used in aggregator).}

\item{deltaD}{numeric, acceptance band for aggregator methods.}
}
\value{
If \code{parameters=FALSE}, a numeric vector of length \code{n} (the random sample).
If \code{parameters=TRUE}, a list with elements:
\item{muX}{the input meanValue}
\item{sdValue}{the fitted standard deviation from partial quantiles (\code{sdQ()})}
\item{X}{the random sample of length \code{n}}
}
\description{
Generate a random sample from a partial-quantile approach, optionally returning fitted parameters.
}
