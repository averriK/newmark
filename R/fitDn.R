# nolint start
#' @title Seismic Newmark Displacement Library
#' @description Core routine for computing permanent Newmark displacements (Dn)
#'   with empirical sliding‑block models and Monte‑Carlo sampling.
#'
#' @details
#' \code{fitDn()} computes weighted quantiles of Newmark displacement from multiple
#' sliding‑block models, optionally sampling hazard (Sa) and/or displacement
#' model uncertainty.
#' (Only the logic errors in “uncertainty” parsing and the quantile call are
#' fixed here.)
#'
#' @param uhs   data.table – uniform‑hazard spectrum (Tn, Sa, p).
#' @param ky    numeric vector of yield accelerations (g).
#' @param Ts    numeric scalar – fundamental period of the sliding mass (s).
#' @param Mw    numeric scalar – scenario moment magnitude (default 6.5).
#' @param NS    integer ≥ 1 – Monte‑Carlo samples per model (default 30).
#' @param models character vector – model identifiers.
#' @param score  numeric vector of same length as models; converted to weights.
#' @param BM_model character, either "crustal"/"shallow" or "interface"/"subduction".
#' @param uncertainty character: any combination of "none", "Sa", "Dn", "both".
#'
#' @return data.table with columns:
#'   \itemize{
#'     \item p   Probability label, e.g. "0.16", "0.5", "mean"
#'     \item Dn  Weighted displacement quantile (cm)
#'   }
#' @import data.table
#' @importFrom Hmisc wtd.quantile
#' @importFrom stats rnorm runif approx
#' @export
fitDn <- function(
    uhs,
    ky,
    Ts,
    Mw = 6.5,
    NS = 30,
    models = c("YG91", "AM88", "JB07", "BT07", "SR08", "BM17", "BM19"),
    score = c(1, 2, 2, 3, 3, 4, 4),
    BM_model = "crustal",
    uncertainty = "none") {
    ## ------------------------------------------------------------------ ##
    ## 1)  **PATCH** – interpret `uncertainty` in a fully case‑insensitive
    ##     and vector‑friendly way.  Nothing else in the function changes.
    ## ------------------------------------------------------------------ ##
    unc_set <- unique(tolower(uncertainty))
    hasSa <- any(unc_set == "sa")
    hasDn <- any(unc_set == "dn")
    unc_mode <- if (hasSa && hasDn) {
        "both"
    } else if (hasSa) {
        "sa"
    } else if (hasDn) {
        "dn"
    } else {
        "none"
    }

    ## ---------- 2. normalise model weights ---------------------------------
    weights <- data.table::data.table(ID = models, weight = score / sum(score))

    ## ---------- 3. pre‑process UHS -----------------------------------------
    UHS <- data.table::copy(uhs)
    UHS[Tn == 0, Tn := 0.01, by = .(p)] # avoid log(0)

    ## ---------- Helper for mean hazard -------------------------------------
    getMeanSa <- function(UHS, Td) {
        tmpMean <- UHS[Tn == Td & p == "mean"]
        if (nrow(tmpMean) == 1L) {
            val <- tmpMean$Sa
        } else if (nrow(tmpMean) == 0) {
            tmp <- UHS[p == "mean"]
            if (nrow(tmp) == 0) stop("No 'mean' row found for hazard.")
            val <- approx(
                x    = log(tmp$Tn),
                y    = log(tmp$Sa),
                xout = log(Td),
                rule = 2
            )$y |> exp()
        } else {
            stop("Multiple 'mean' rows for Tn = ", Td)
        }
        rep(val, NS)
    }

    ## ---------- 4. sample Sa at required periods ---------------------------
    if (unc_mode %in% c("sa", "both")) {
        PGA <- sampleSa(UHS, Td = 0.01, n = NS)$Sa
        Sa_10_Ts <- sampleSa(UHS, Td = 1.0 * Ts, n = NS)$Sa
        Sa_13_Ts <- sampleSa(UHS, Td = 1.3 * Ts, n = NS)$Sa
        Sa_15_Ts <- sampleSa(UHS, Td = 1.5 * Ts, n = NS)$Sa
    } else {
        PGA <- getMeanSa(UHS, 0.01)
        Sa_10_Ts <- getMeanSa(UHS, 1.0 * Ts)
        Sa_13_Ts <- getMeanSa(UHS, 1.3 * Ts)
        Sa_15_Ts <- getMeanSa(UHS, 1.5 * Ts)
    }

    PGV <- (PGA^1.0529) * exp(0.1241) * 100 # [cm/s]
    AI <- (PGA^1.9228) * exp(2.6109) # [m/s]

    ## ---------- 5. Monte‑Carlo sampling per model --------------------------
    # ---------------------------------------------------------------------------
    # Helper: draw Monte‑Carlo samples for one empirical Dn model
    # ---------------------------------------------------------------------------
    getDnUncertainty <- function(DnFun, ..., n = 1) {
        DnModel <- DnFun(...) # data.table(muLnD, sdLnD, ID)
        stopifnot(all(c("muLnD", "sdLnD", "ID") %in% names(DnModel)))

        ## --- PATCH: never call rnorm() with NA / non‑finite arguments ----------


        if (unc_mode %in% c("dn", "both")) {
            out <- DnModel[
                , .(ID, LnD = safe_rnorm(n, muLnD, sdLnD)), # <- patched
                by = .I
            ][
                , .(sample = .I, ID, Dn = exp(LnD))
            ]
        } else { # deterministic branch
            out <- DnModel[
                , .(sample = seq_len(n), ID, Dn = exp(muLnD)),
                by = .I
            ]
        }
        out
    }


    DnTable <- data.table::data.table()
    add <- function(dt, new) data.table::rbindlist(list(dt, new), use.names = TRUE, fill = TRUE)

    if ("AM88" %in% models) {
        DnTable <- add(DnTable, getDnUncertainty(Dn_AM88, PGA = PGA, ky = ky, n = NS))
    }
    if ("YG91" %in% models) {
        DnTable <- add(DnTable, getDnUncertainty(Dn_YG91, PGA = PGA, ky = ky, n = NS))
    }
    if ("JB07" %in% models) {
        DnTable <- add(DnTable, getDnUncertainty(Dn_JB07, PGA = PGA, AI = AI, ky = ky, n = NS))
    }
    if ("SR08" %in% models) {
        DnTable <- add(DnTable, getDnUncertainty(Dn_SR08, PGA = PGA, AI = AI, ky = ky, n = NS))
    }
    if ("BT07" %in% models) {
        DnTable <- add(DnTable, getDnUncertainty(Dn_BT07, Ts = Ts, Sa = Sa_10_Ts, Mw = Mw, ky = ky, n = NS))
    }
    if (tolower(BM_model) %in% c("interface", "subduction") && "BM17" %in% models) {
        DnTable <- add(DnTable, getDnUncertainty(Dn_BM17, Ts = Ts, Sa = Sa_15_Ts, Mw = Mw, ky = ky, n = NS))
    }
    if (tolower(BM_model) %in% c("shallow", "crustal") && "BM19" %in% models) {
        DnTable <- add(DnTable, getDnUncertainty(Dn_BM19,
            Ts = Ts, Sa = Sa_13_Ts,
            PGA = PGA, PGV = PGV, Mw = Mw,
            ky = ky, n = NS
        ))
    }

    ## ---------- 6. Weighted quantiles (FIX #2 — type must be integer 1‑9) ---
    AUX <- weights[DnTable, on = "ID"]
    AUX[, units := "cm"]

    probs <- UHS[p != "mean", unique(as.numeric(p))]

    Dn_q <- AUX[, .(
        p = probs,
        Dn = wtd.quantile(
            x       = Dn,
            weights = weight,
            probs   = probs,
            # type    = 7, # ← numeric, not string: obsoleted in Hmisc 5.0.0
            type    = "quantile", # ← misc 5.1
            na.rm   = TRUE
        )
    )]

    Dn_mean <- data.table::data.table(p = "mean", Dn = mean(AUX$Dn))

    data.table::rbindlist(list(Dn_q, Dn_mean))
}

safe_rnorm <- function(n, m, s) {
    if (is.finite(m) && is.finite(s) && s > 0) {
        stats::rnorm(n, m, s)
    } else { # outside model validity ⇒ return NA
        rep(NA_real_, n)
    }
}
