% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitDn.R
\name{fitDn}
\alias{fitDn}
\title{Seismic Newmark Displacement Library}
\usage{
fitDn(
  uhs,
  ky,
  Ts,
  Mw = 6.5,
  NS = 30,
  models = c("YG91", "AM88", "JB07", "BT07", "SR08", "BM17", "BM19"),
  score = c(1, 2, 2, 3, 3, 4, 4),
  BM_model = "crustal",
  uncertainty = "none"
)
}
\arguments{
\item{uhs}{data.table – uniform-hazard spectrum (Tn, Sa, p).}

\item{ky}{numeric vector of yield accelerations (g).}

\item{Ts}{numeric scalar – fundamental period of the sliding mass (s).}

\item{Mw}{numeric scalar – scenario moment magnitude (default 6.5).}

\item{NS}{integer ≥ 1 – Monte-Carlo samples per model (default 30).}

\item{models}{character vector – model identifiers.}

\item{score}{numeric vector of same length as models; converted to weights.}

\item{BM_model}{character, either "crustal"/"shallow" or "interface"/"subduction".}

\item{uncertainty}{character: any combination of "none", "Sa", "Dn", "both".
For example, user can enter \code{c("sa","dn")} or just "sA", etc.
\itemize{
\item "none": no hazard fractile sampling, no displacement model random draws
\item "Sa": hazard fractile sampling only
\item "Dn": displacement model random draws only
\item "both": hazard + displacement model random draws
}}
}
\value{
data.table with columns:
\itemize{
\item p (probability label, e.g. "0.16", "0.5", "mean")
\item Dn (cm) Weighted quantiles from all models
}
}
\description{
Core routines for computing permanent Newmark displacements (Dn)
with empirical sliding-block models and Monte-Carlo sampling.
}
\details{
\code{fitDn()} computes weighted quantiles of Newmark displacement from multiple
sliding-block models, optionally sampling hazard (Sa) and/or displacement
model uncertainty.
}
