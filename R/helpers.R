# R/helpers.R

sample_Sa <- function(UHS, Td, n) {
    if (Td %in% UHS$Tn) {
        # Select
        SaTable <- UHS[Tn == Td, .(Sa, p)]
    }
    if (!(Td %in% UHS$Tn)) {
        # Interpolate
        SaTable <- UHS[, .(Sa = stats::approx(
            rule = 2,
            x = log(Tn), y = log(Sa), xout = log(Td)
        )$y |> exp()), by = .(p)]
    }
    muLnSa <- log(SaTable[p == "mean"]$Sa)
    p <- as.numeric(SaTable[p != "mean"]$p)
    q <- log(SaTable[p != "mean"]$Sa)
    sdLnSa <- sdQ(meanValue = muLnSa, p = p, q = q)
    Sa <- rnormQ(meanValue = muLnSa, p = p, q = q, n = n) |> exp()
    return(list(SaTable = SaTable, muLnSa = muLnSa, sdLnSa = sdLnSa, Sa = Sa))
}



