% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitSaF.R
\name{fitSaF}
\alias{fitSaF}
\title{Site-Factor Amplification with Uncertainty}
\usage{
fitSaF(
  uhs,
  vs30,
  vref = 760,
  NS = 30,
  models = c("ST17"),
  score = rep(1, length(models)),
  uncertainty = c("none", "Sa", "F", "both")
)
}
\arguments{
\item{uhs}{data.table with columns (Tn, p, Sa). Must include p="mean" plus
0<p<1 rows if you do hazard sampling. Must also have Tn=0 or Tn=0.01 for PGA
if your site-factor model depends on that.}

\item{vs30}{numeric scalar – target site velocity.}

\item{vref}{numeric – reference velocity (760 by default).}

\item{NS}{integer – number of Monte Carlo draws (default 30).}

\item{models}{character – site-factor models, default "ST17".}

\item{score}{numeric – weighting for the model(s).}

\item{uncertainty}{character vector – any combination of c("Sa","F"):
\itemize{
\item "none" => no hazard fractiles, no site-factor random
\item "Sa"   => hazard fractiles only
\item "F"    => site-factor random only
\item "both" => hazard + site-factor random
}}
}
\value{
data.table with columns:
\itemize{
\item Tn  => oscillator period
\item p   => probabilities from the original p set (excluding p="mean")
\item SaF => weighted quantiles of final \eqn{Sa * F}
}
plus a "mean" row for p="mean".
}
\description{
Similar to \code{fitDn()} but for site-amplified spectral acceleration (SaF).
For each oscillator period Tn in the input table, we sample hazard (if "Sa" in uncertainty),
compute or sample site factor (if "F" in uncertainty), then produce draws of \eqn{SaF = Sa * F}.
Finally, we aggregate them via weighted quantiles in the p dimension.
}
