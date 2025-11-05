# nolint start
#' Check that POE(p) does not decrease (allows ties within tol.)
#' @param AEPTable data.table with columns lon,lat,depth,Tn,Sa,p,POE
#' @param tolAbs   absolute tolerance (default 1e‑8)
#' @param tolRel   relative tolerance with respect to POE (default 1e‑4)
#' @param strict   TRUE -> stop if fails; FALSE -> only message.
#'
#' @return TRUE invisible; or data.table with problematic curves when strict = FALSE.
#' @export
checkAEP <- function(AEPTable,
                     tolAbs = 1e-8,
                     tolRel = 1e-4,
                     strict = TRUE) {
    stopifnot(inherits(AEPTable, "data.table"))
    DT <- AEPTable[p != "mean"] # quantiles
    if (!nrow(DT)) {
        return(invisible(TRUE))
    }

    bad <- DT[,
        {
            poe <- POE[order(p)]
            DELTA <- diff(poe)
            ok <- DELTA >= -pmax(tolAbs, tolRel * poe[-length(poe)])
            if (any(!ok)) .SD
        },
        by = .(lon, lat, depth, Tn, Sa)
    ]

    if (nrow(bad) && strict) {
        stop(
            sprintf(
                "checkAEP(): %d non-monotonic curves.",
                uniqueN(bad, by = .(lon, lat, depth, Tn, Sa))
            ),
            call. = FALSE
        )
    }
    if (nrow(bad) && !strict) {
        message(sprintf(
            "checkAEP(): %d non-monotonic curves (only warning).",
            uniqueN(bad, by = .(lon, lat, depth, Tn, Sa))
        ))
        return(bad)
    }
    invisible(TRUE)
}
