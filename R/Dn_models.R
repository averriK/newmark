# nolint start
#  Displacement models (empirical equations) --------------------------------
#' Ambraseys & Menu (1988) rigid sliding-block model
#' @param PGA Peak ground acceleration (g).
#' @param ky  Yield acceleration (g).
#' @return data.table(muLnD, sdLnD, ID = "AM88")
#' @export
Dn_AM88 <- function(PGA, ky) {
    r <- ky / PGA
    r <- pmin(r, 0.9999) # keep (1 - r) > 0
    muLog10D <- 0.90 + log10((1 - r)^2.53 * r^-1.09)
    sdLog10D <- 0.30

    data.table(
        muLnD = muLog10D * log(10),
        sdLnD = sdLog10D * log(10),
        ID = "AM88"
    )
}

#' Yegian et al. (1991) rigid sliding-block model (patched to handle ky ≥ PGA)
#'
#' @inheritParams Dn_AM88
#' @return data.table with columns:
#'   \describe{
#'     \item{muLnD}{log-mean Newmark displacement (natural log [cm])}
#'     \item{sdLnD}{log-stddev Newmark displacement}
#'     \item{ID}{model identifier, "YG91"}
#'   }
#' @export
Dn_YG91 <- function(PGA, ky) {
    # ratio of yield to demand
    r <- ky / PGA

    # compute polynomial for r < 1
    muLog10D <- 0.22 - 10.12 * r + 16.38 * r^2 - 11.48 * r^3

    sdLog10D <- 0.45

    # convert to natural-log scale and return
    data.table(
        muLnD = muLog10D * log(10),
        sdLnD = sdLog10D * log(10),
        ID = "YG91"
    )
}


#' Jibson (2007) empirical model
#' @inheritParams Dn_AM88
#' @param AI Arias Intensity (m/s); if NULL, back-calculated from PGA.
#' @return data.table(muLnD, sdLnD, ID = "JB07")
#' @export
Dn_JB07 <- function(PGA, ky, AI = NULL) {
    if (is.null(AI)) AI <- (PGA^1.9228) * exp(2.6109)
    r <- ky / PGA
    muLog10D <- 0.561 * log10(AI) - 3.833 * log10(r) - 1.474
    sdLog10D <- 0.616

    data.table(
        muLnD = muLog10D * log(10),
        sdLnD = sdLog10D * log(10),
        ID = "JB07"
    )
}

#' Saygili & Rathje (2008) model
#' @inheritParams Dn_JB07
#' @return data.table(muLnD, sdLnD, ID = "SR08")
#' @export
Dn_SR08 <- function(PGA, ky, AI = NULL) {
    if (is.null(AI)) AI <- (PGA^1.9228) * exp(2.6109)
    r <- ky / PGA
    muLnD <- 2.39 - 5.24 * r - 18.78 * r^2 + 42.01 * r^3 - 29.15 * r^4 -
        1.56 * log(PGA) + 1.38 * log(AI)
    sdLnD <- 0.46 + 0.56 * r

    data.table(
        muLnD = muLnD,
        sdLnD = sdLnD,
        ID = "SR08"
    )
}

#' Bray & Travasarou (2007) flexible sliding-block model
#' @param ky Yield acceleration (g).
#' @param Sa Spectral acceleration at 1.5 * Ts (g).
#' @param Ts Fundamental period (s).
#' @param Mw Moment magnitude.
#' @return data.table(muLnD, sdLnD, ID = "BT07")
#' @export
Dn_BT07 <- function(ky, Sa, Ts, Mw = 6.5) {
    lnSa <- log(Sa)
    lnky <- log(ky)

    muLnD <- -1.10 - 2.83 * lnky - 0.333 * lnky^2 +
        0.566 * lnky * lnSa + 3.04 * lnSa -
        0.244 * lnSa^2 + 0.278 * (Mw - 7)

    sdLnD <- 0.66

    data.table(
        muLnD = muLnD,
        sdLnD = sdLnD,
        ID = "BT07"
    )
}

#' Bray & Macedo (2017) – subduction/interface events
#' @param ky Yield acceleration (g).
#' @param Sa Spectral acceleration at 1.5 * Ts (g).
#' @param Ts Fundamental period (s).
#' @param Mw Moment magnitude (default 7.5).
#' @return data.table(muLnD, sdLnD, ID = "BM17")
#' @export
Dn_BM17 <- function(ky, Sa, Ts, Mw = 7.5) {
    lnSa <- log(Sa)

    if (Ts < 0.1) {
        c1 <- -5.864
        c2 <- -9.421
        c3 <- 0.0
    } else {
        c1 <- -6.896
        c2 <- 3.081
        c3 <- -0.803
    }

    a0 <- c1 + 0.550 * Mw + c2 * Ts + c3 * Ts^2 -
        3.353 * log(ky) - 0.390 * log(ky)^2
    a1 <- 3.060 + 0.538 * log(ky)
    a2 <- -0.225

    muLnD <- a0 + a1 * lnSa + a2 * lnSa^2
    sdLnD <- 0.73

    data.table(
        muLnD = muLnD,
        sdLnD = sdLnD,
        ID = "BM17"
    )
}

#' Bray & Macedo (2019) – shallow-crustal events
#' @param ky  Yield acceleration (g).
#' @param Ts  Fundamental period (s).
#' @param Sa  Spectral acceleration at 1.3 * Ts (g).
#' @param PGA Peak ground acceleration (g).
#' @param PGV Peak ground velocity (cm/s).
#' @param Mw  Moment magnitude (default 6.5).
#' @return data.table(muLnD, sdLnD, ID = "BM19")
#' @export
Dn_BM19 <- function(ky, Ts, Sa, PGA, PGV, Mw = 6.5) {
    n <- length(Sa)
    if (!(length(PGA) == n && length(PGV) == n)) {
        stop("Sa, PGA, and PGV must have the same length")
    }

    lnSa <- log(Sa)
    I <- PGV <= 115 # TRUE = ordinary, FALSE = pulse

    c1 <- ifelse(Ts < 0.1, -4.551, -5.894)
    c2 <- ifelse(Ts < 0.1, -9.690, 3.152)
    c3 <- ifelse(Ts < 0.1, 0.000, -0.910)

    lnPGV_term <- ifelse(I, 0, 1)
    PGV_shift <- ifelse(I, 0, -log(115))

    a0 <- c1 + 0.607 * Mw + c2 * Ts + c3 * Ts^2 -
        2.491 * log(ky) - 0.245 * log(ky)^2 +
        lnPGV_term * log(PGV) + PGV_shift

    a1 <- 2.703 + 0.344 * log(ky)
    a2 <- -0.089

    muLnD <- a0 + a1 * lnSa + a2 * lnSa^2
    sdLnD <- 0.74

    data.table(
        muLnD = muLnD,
        sdLnD = sdLnD,
        ID = "BM19"
    )
}
