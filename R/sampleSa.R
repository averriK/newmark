# nolint start
#' Sample spectral acceleration consistent with UHS quantiles
#' @param UHS data.table with Tn, Sa, p (prob labels + "mean").
#' @param Td  numeric scalar – target period (s).
#' @param n   integer – number of samples.
#' @return list(SaTable, muLnSa, sdLnSa, Sa)
#' @import data.table
#' @importFrom sdQ rnormQ sdQ
#' @export
sampleSa <- function(UHS, Td, n) {
    if (Td %in% UHS$Tn) {
        SaTable <- UHS[Tn == Td, .(Sa, p)]
    } else {
        SaTable <- UHS[
            , .(Sa = stats::approx(
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



# ---------------------------------------------------------------------------

# nolint end
