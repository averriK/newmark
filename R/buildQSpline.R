# nolint start
#' Build a monotone quantile spline Q(u) for ln(Sa)
#' @importFrom stats splinefun pnorm
#' @importFrom mvtnorm rmvnorm
#' @param SaTable data.table with probability column \code{p} (0-1 or "mean")
#'                and spectral-acceleration column \code{Sa}.
#' @return A closure \code{Q(u)} returning \code{ln(Sa)} for any \code{u in[0,1]}.
#' @export
buildQSpline <- function(SaTable) {
    requiredCols <- c("p", "Sa")
    if (!all(requiredCols %in% names(SaTable))) {
        stop("buildQSpline(): SaTable must contain columns 'p' and 'Sa'.") # nolint
    }

    DT <- data.table::as.data.table(SaTable)

    # discard 'mean' row
    DT <- DT[p != "mean"]

    if (nrow(DT) < 3L) {
        stop("buildQSpline(): 3 or more quantiles are required.") # nolint
    }

    extraCols <- setdiff(names(DT), c("p", "Sa"))
    if (length(extraCols)) {
        bad <- extraCols[vapply(
            DT[, extraCols, with = FALSE],
            function(x) data.table::uniqueN(x) > 1L,
            logical(1)
        )]
        if (length(bad)) {
            stop(
                "buildQSpline(): non-unique columns: ", # nolint
                paste(bad, collapse = ", ")
            )
        }
    }

    pVec <- sort(as.numeric(DT$p))
    lnSaVec <- log(DT[order(as.numeric(p)), Sa])

    lnSaMono <- cummax(lnSaVec) # avoid descents
    if (any(diff(lnSaMono) <= 0)) { # break ties exactly
        lnSaMono <- lnSaMono + seq_along(lnSaMono) * 1e-8
    }

    splineFun <- stats::splinefun(
        x = pVec,
        y = lnSaMono,
        method = "hyman"
    )

    Q <- function(u) {
        u <- pmin(pmax(u, 0), 1) # clamp to [0,1]
        splineFun(u)
    }

    Q
}
