# ---------------------------------------------------------------------------
#' Log‑normal sampling of spectral acceleration consistent with UHS quantiles
#'
#' @description
#' Generates a Monte‑Carlo realisation of spectral acceleration (Sa) at a target
#' period \code{Td}, consistent with the quantile structure of a uniform‑hazard
#' spectrum (UHS).  Dispersion is inferred via the helper [sdQ()], and samples
#' are drawn with [rnormQ()].
#'
#' @param UHS \link[data.table]{data.table} with columns \code{Tn}, \code{Sa},
#'   and \code{p}.  Includes a special row where \code{p == "mean"}.  See
#'   [fitDn()].
#' @param Td Numeric scalar. Target oscillator period (s).
#' @param n Integer \(\ge 1\). Number of random samples to draw.
#'
#' @return List with elements
#' \describe{
#'   \item{SaTable}{\link[data.table]{data.table} of Sa quantiles at \code{Td}.}
#'   \item{muLnSa}{Mean of \eqn{\ln Sa}.}
#'   \item{sdLnSa}{Standard deviation of \eqn{\ln Sa}.}
#'   \item{Sa}{Numeric vector of length \code{n} with sampled Sa values (g).}
#' }
#' @export
sample_Sa <- function(UHS, Td, n) {
    if (!inherits(UHS, "data.table")) UHS <- data.table::as.data.table(UHS)
    if (Td %in% UHS$Tn) {
        SaTable <- UHS[Tn == Td, .(Sa, p)]
    } else {
        SaTable <- UHS[, .(Sa = stats::approx(
            x = log(Tn),
            y = log(Sa),
            xout = log(Td),
            rule = 2
        )$y |> exp()),
        by = .(p)
        ]
    }
    muLnSa <- log(SaTable[p == "mean", Sa])
    p_vec <- as.numeric(SaTable[p != "mean", p])
    q_vec <- log(SaTable[p != "mean", Sa])
    sdLnSa <- sdQ(meanValue = muLnSa, p = p_vec, q = q_vec)
    Sa <- rnormQ(meanValue = muLnSa, p = p_vec, q = q_vec, n = n) |> exp()
    list(SaTable = SaTable, muLnSa = muLnSa, sdLnSa = sdLnSa, Sa = Sa)
}
