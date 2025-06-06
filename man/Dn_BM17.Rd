% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Dn_models.R
\name{Dn_BM17}
\alias{Dn_BM17}
\title{Bray & Macedo (2017) – subduction/interface events}
\usage{
Dn_BM17(ky, Sa, Ts, Mw = 7.5)
}
\arguments{
\item{ky}{Yield acceleration (g).}

\item{Sa}{Spectral acceleration at 1.5 * Ts (g).}

\item{Ts}{Fundamental period (s).}

\item{Mw}{Moment magnitude (default 7.5).}
}
\value{
data.table(muLnD, sdLnD, ID = "BM17")
}
\description{
Bray & Macedo (2017) – subduction/interface events
}
