# ============================================================
#  sdSaQ()    –  constant ln-σ* (median of period MAD values)
# ============================================================
#' Robust period-independent ln-σ for an Sa spectrum
#'
#' @param UHS  data.table with Tn, Sa, p  (must include row p == "mean").
#' @return numeric scalar  —  σ*  (ln units).
#' @keywords internal
sdSaQ <- function(UHS) {
    stopifnot(all(c("Tn", "Sa", "p") %in% names(UHS)))

    stats::median(
        UHS[p != "mean",
            {
                z <- log(Sa)
                1.4826 * stats::median(abs(z - stats::median(z)))
            },
            by = .(Tn)
        ]$V1
    )
}
