# nolint start
#' Check monotonicity and duplicates in UHSTable
#'
#' Ensures that, for every combination
#' \code{TR, Vs30, Vref, ID, Tn}, there is exactly
#' one row per probability \code{p} and that
#' \code{Sa(p)} is non‑decreasing.
#'
#' @param UHSTable data.table with columns
#'   \code{TR, Vs30, Vref, ID, Tn, p, Sa}.
#' @param epsMonot numeric tolerance (default 1e-6).
#'
#' @return invisibly \code{TRUE}; stops on error.
#' @export
checkUHS <- function(UHSTable, epsMonot = 1e-6) {
    stopifnot(inherits(UHSTable, "data.table"))

    dup <- UHSTable[p != "mean",
        .N,
        by = .(TR, Vs30, Vref, ID, Tn, p)
    ][N > 1L]
    if (nrow(dup)) {
        stop(
            sprintf(
                "checkUHS(): %d duplicados de cuantil.",
                nrow(dup)
            ),
            call. = FALSE
        )
    }

    bad <- UHSTable[p != "mean",
        {
            s <- Sa[order(p)]
            if (any(diff(s) < -epsMonot)) .SD
        },
        by = .(TR, Vs30, Vref, ID, Tn)
    ]

    if (nrow(bad)) {
        stop(
            sprintf(
                "checkUHS(): ERROR",
                uniqueN(bad, by = .(TR, Vs30, Vref, ID, Tn))
            ),
            call. = FALSE
        )
    }

    invisible(TRUE)
}
