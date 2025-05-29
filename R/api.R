# nolint start
# ---------------------------------------------------------------------------
#' @title Seismic Newmark Displacement Library
#' @description Core routines for computing permanent Newmark displacements (Dn)
#'   with empirical sliding-block models and Monte-Carlo sampling.
#' @import data.table
#' @importFrom Hmisc wtd.quantile
#' @keywords internal
NULL
# ---------------------------------------------------------------------------

#' Weighted Newmark displacement quantiles from multiple models
#'
#' @param UHS   data.table – uniform-hazard spectrum (cols: Tn, Sa, p).
#' @param ky    numeric vector of yield accelerations (g).
#' @param Ts    numeric scalar – fundamental period of the sliding mass (s).
#' @param Mw    numeric scalar – scenario moment magnitude (default 6.5).
#' @param NS    integer ≥ 1 – Monte-Carlo sample size per model (default 30).
#' @param models character vector – model identifiers.
#' @param score  numeric vector of same length as models; converted to weights.
#' @param BM_model "crustal"/"shallow" or "interface"/"subduction".
#' @return data.table with probability label p and weighted displacement Dn (cm).
#' @export
fitDn <- function(
    UHS,
    ky,
    Ts,
    Mw = 6.5,
    NS = 30,
    models = c("YG91", "AM88", "JB07", "BT07", "SR08", "BM17", "BM19"),
    score = c(1, 2, 2, 3, 3, 4, 4),
    BM_model = "crustal") {
    weights <- data.table(ID = models, weight = score / sum(score))

    ## --- pre-process UHS ----------------------------------------------------
    UHS[Tn == 0, Tn := 0.01, by = .(p)] # avoid log(0)

    ## --- Sa samples at required periods ------------------------------------
    PGA <- sample_Sa(UHS, Td = 0.01, n = NS)$Sa
    PGV <- (PGA^1.0529) * exp(0.1241) * 100 # [cm/s]
    AI <- (PGA^1.9228) * exp(2.6109) # [m/s]
    Sa_10_Ts <- sample_Sa(UHS, Td = 1.0 * Ts, n = NS)$Sa
    Sa_13_Ts <- sample_Sa(UHS, Td = 1.3 * Ts, n = NS)$Sa
    Sa_15_Ts <- sample_Sa(UHS, Td = 1.5 * Ts, n = NS)$Sa

    ## --- Monte-Carlo sampling per model ------------------------------------
    DnTable <- data.table()
    append <- function(dt, new) rbindlist(list(dt, new), use.names = TRUE, fill = TRUE)

    if ("AM88" %in% models) {
        DnTable <- append(DnTable, Dn_model(Dn_AM88, PGA = PGA, ky = ky, n = NS))
    }
    if ("YG91" %in% models) {
        DnTable <- append(DnTable, Dn_model(Dn_YG91, PGA = PGA, ky = ky, n = NS))
    }
    if ("JB07" %in% models) {
        DnTable <- append(DnTable, Dn_model(Dn_JB07, PGA = PGA, AI = AI, ky = ky, n = NS))
    }
    if ("SR08" %in% models) {
        DnTable <- append(DnTable, Dn_model(Dn_SR08, PGA = PGA, AI = AI, ky = ky, n = NS))
    }
    if ("BT07" %in% models) {
        DnTable <- append(
            DnTable,
            Dn_model(Dn_BT07, Ts = Ts, Sa = Sa_10_Ts, Mw = Mw, ky = ky, n = NS)
        )
    }
    if (tolower(BM_model) %in% c("interface", "subduction") && "BM17" %in% models) {
        DnTable <- append(
            DnTable,
            Dn_model(Dn_BM17, Ts = Ts, Sa = Sa_15_Ts, Mw = Mw, ky = ky, n = NS)
        )
    }
    if (tolower(BM_model) %in% c("shallow", "crustal") && "BM19" %in% models) {
        DnTable <- append(
            DnTable,
            Dn_model(
                Dn_BM19,
                Ts  = Ts,
                Sa  = Sa_13_Ts,
                PGA = PGA,
                PGV = PGV,
                Mw  = Mw,
                ky  = ky,
                n   = NS
            )
        )
    }

    ## --- Weighted quantiles -------------------------------------------------
    AUX <- weights[DnTable, on = "ID"]
    probs <- sort(unique(as.numeric(UHS[p != "mean", p])))

    Dn_q <- AUX[
        , .(
            p = probs,
            Dn = Hmisc::wtd.quantile(Dn,
                weights = weight,
                probs = probs,
                type = "quantile",
                na.rm = TRUE
            )
        ),
        by = .(dummy = 1)
    ][, dummy := NULL]

    Dn_mean <- data.table(p = "mean", Dn = mean(AUX$Dn))
    rbindlist(list(Dn_q, Dn_mean))[, units := "cm"]
}
