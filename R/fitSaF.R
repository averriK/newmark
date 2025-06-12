#' @title Site-Factor Amplification with Uncertainty
#'
#' @description Similar to \code{fitDn()} but for site-amplified spectral acceleration (SaF).
#'   For each oscillator period Tn in the input table, we sample hazard (if "Sa" in uncertainty),
#'   compute or sample site factor (if "F" in uncertainty), then produce draws of \eqn{SaF = Sa * F}.
#'   Finally, we aggregate them via weighted quantiles in the p dimension.
#'
#' @param uhs data.table with columns (Tn, p, Sa). Must include p="mean" plus
#'   0<p<1 rows if you do hazard sampling. Must also have Tn=0 or Tn=0.01 for PGA
#'   if your site-factor model depends on that.
#' @param vs30 numeric scalar – target site velocity.
#' @param vref numeric – reference velocity (760 by default).
#' @param NS integer – number of Monte Carlo draws (default 30).
#' @param models character – site-factor models, default "ST17".
#' @param score numeric – weighting for the model(s).
#' @param uncertainty character vector – any combination of c("Sa","F"):
#'   * "none" => no hazard fractiles, no site-factor random
#'   * "Sa"   => hazard fractiles only
#'   * "F"    => site-factor random only
#'   * "both" => hazard + site-factor random
#'
#' @return data.table with columns:
#'   \itemize{
#'     \item Tn  => oscillator period
#'     \item p   => probabilities from the original p set (excluding p="mean")
#'     \item SaF => weighted quantiles of final \eqn{Sa * F}
#'   }
#'   plus a "mean" row for p="mean".
#' @import data.table
#' @importFrom stats rnorm
#' @export
fitSaF <- function(
    uhs,
    vs30,
    vref = 760,
    NS = 30,
    models = c("ST17"),
    score = rep(1, length(models)),
    uncertainty = c("none", "Sa", "F", "both")) {
    ## ------------------------------------------------------------------------
    ## 1) Parse `uncertainty`
    ## ------------------------------------------------------------------------
    unc_vec <- tolower(uncertainty)
    doSa <- any(unc_vec %in% "sa")
    doF <- any(unc_vec %in% "f")
    mode <- "none"
    if (doSa && doF) {
        mode <- "both"
    } else if (doSa && !doF) {
        mode <- "sa"
    } else if (!doSa && doF) {
        mode <- "f"
    }

    ## ------------------------------------------------------------------------
    ## 2) Build model-weight table. Insert "gem" if not present
    ## ------------------------------------------------------------------------
    AUX <- data.table(ID = models, weight = score)
    if (!("gem" %in% AUX$ID)) {
        AUX <- rbind(AUX, data.table(ID = "gem", weight = 1))
    }
    AUX[, weight := weight / sum(weight)]
    wTable <- data.table::copy(AUX) # store final weight table

    ## ------------------------------------------------------------------------
    ## 3) Preprocess UHS: Tn=0 => Tn=0.01
    ## ------------------------------------------------------------------------
    UHS <- data.table::copy(uhs)
    UHS[Tn == 0, Tn := 0.01, by = .(p)]

    # Helper to sample hazard
    getHazardDraws <- function(UHS, period, n = 1) {
        if (mode %in% c("sa", "both")) {
            # aggregator from partial-quantile or piecewise
            out <- sampleSa(UHS, Td = period, n = n)$Sa
        } else {
            # replicate the mean
            TEMP <- UHS[Tn == period & p == "mean"]
            if (nrow(TEMP) != 1) {
                stop("No unique mean row at Tn=", period)
            }
            out <- rep(TEMP$Sa, times = n)
        }
        out
    }

    ## ------------------------------------------------------------------------
    ## 4) We'll produce a final aggregator table
    ## ------------------------------------------------------------------------
    RES <- data.table() # final aggregator across Tn

    # gather Tn>0.01
    TnList <- sort(unique(UHS[Tn > 0.009, Tn]))

    # sample pga once at Tn=0.01
    pga <- getHazardDraws(UHS, period = 0.01, n = NS)

    for (thisTn in TnList) {
        # hazard draws for that Tn
        SaRock <- getHazardDraws(UHS, period = thisTn, n = NS)

        # Build a single table named AUX that stores sample, ID, SaRock, SaF
        if (vs30 == vref) {
            # no amplification => ID="gem"
            AUX <- data.table(
                sample = 1:NS,
                ID = "gem",
                SaRock = SaRock,
                SaF = SaRock # same as rock
            )
        } else {
            # vs30 != vref => aggregator with site-factor models
            AUX <- data.table()
            for (modID in models) {
                TEMP <- data.table(
                    sample = 1:NS,
                    ID     = modID,
                    SaRock = SaRock,
                    pga    = pga
                )

                # row-by-row aggregator
                # if "F" => random => else deterministic
                TEMP[, SaF := {
                    # e.g. call SaF_ST17
                    tmp <- SaF_ST17(
                        Sa = SaRock,
                        pga = pga,
                        Tn = thisTn,
                        vs30 = vs30,
                        vref = vref
                    )
                    if (mode %in% c("f", "both")) {
                        # random
                        LnSaF <- rnorm(1, tmp$muLnSaF, tmp$sdLnSaF)
                        exp(LnSaF)
                    } else {
                        # deterministic
                        exp(tmp$muLnSaF)
                    }
                }, by = .(sample)]

                AUX <- rbind(AUX, TEMP[, .(sample, ID, SaRock, SaF)])
            }
        }

        # aggregator => join with wTable => get weight
        AUX <- wTable[AUX, on = "ID"]

        # define p
        PROBS <- UHS[p != "mean", unique(as.numeric(p))]

        # Weighted quantiles => SaRock
        # (some calls to Hmisc)
        SaRock_q <- AUX[, .(
            p = PROBS,
            SaRock = Hmisc::wtd.quantile(
                x = SaRock,
                weights = weight,
                probs = PROBS,
                # type    = 7, # ← numeric, not string: obsoleted in Hmisc 5.0.0
                type = "quantile", # ← Hmisc 5.1
                na.rm = TRUE
            )
        )]
        # Weighted quantiles => SaF
        SaF_q <- AUX[, .(
            p = PROBS,
            SaF = Hmisc::wtd.quantile(
                x = SaF,
                weights = weight,
                probs = PROBS,
                # type    = 7, # ← numeric, not string: obsoleted in Hmisc 5.0.0
                type = "quantile", # ← Hmisc 5.1
                na.rm = TRUE
            )
        )]

        # add mean row
        mean_rock <- data.table(p = "mean", SaRock = mean(AUX$SaRock))
        mean_site <- data.table(p = "mean", SaF = mean(AUX$SaF))

        # combine side-by-side => cbind
        OUTQ <- cbind(
            SaRock_q[, .(p, SaRock)],
            SaF = SaF_q$SaF
        )
        # add mean row
        MeanRow <- data.table(
            p = "mean",
            SaRock = mean_rock$SaRock,
            SaF = mean_site$SaF
        )
        OUTQ <- rbind(OUTQ, MeanRow)
        OUTQ[, Tn := thisTn]

        RES <- rbind(RES, OUTQ)
    }

    data.table::setcolorder(RES, c("Tn", "p", "SaRock", "SaF"))
    RES[]
}
