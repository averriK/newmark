#  Displacement models (empirical equations) --------------------------------

#' Ambraseys & Menu (1988) rigid sliding-block model
#' @param PGA Peak ground acceleration (g).
#' @param ky  Yield acceleration (g).
#' @return data.table(muLnD, sdLnD, ID = "AM88")
#' @export
Dn_AM88 <- function(PGA, ky) {
  r <- ky / PGA
  muLog10D <- 0.90 + log10((1 - r)^2.53 * r^-1.09)
  sdLog10D <- 0.30
  
  data.table(
    muLnD = muLog10D * log(10),
    sdLnD = sdLog10D * log(10),
    ID = "AM88"
  )
}

#' Yegian et al. (1991) rigid sliding-block model
#' @inheritParams Dn_AM88
#' @return data.table(muLnD, sdLnD, ID = "YG91")
#' @export
Dn_YG91 <- function(PGA, ky) {
  r <- ky / PGA
  muLog10D <- 0.22 - 10.12 * r + 16.38 * r^2 - 11.48 * r^3
  sdLog10D <- 0.45
  
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

# ============================================================
#  SaF model: SaF_ST17()
#  (analogous to Dn_AM88, Dn_YG91, …)
# ============================================================
# nolint start
#' Site–factor model SaF_ST17
#'
#' Computes median (muLnSaF) and log-standard deviation (sdLnSaF) of
#' \deqn{SaF = Sa(T_n)\;×\;F_{ST17}(T_n,\mathrm{PGA},V_{s30})}.
#'
#' @param Sa    numeric vector – sampled spectral acceleration at Tn (g).
#' @param PGA   numeric vector – sampled peak ground acceleration (g).
#' @param Tn    numeric scalar – oscillator period (s).
#' @param vs30  numeric scalar – time-averaged Vs30 (m/s).
#' @param vref  numeric scalar – reference Vs (default 760 m/s).
#'
#' @return data.table(muLnSaF, sdLnSaF, ID = "ST17")
#' @export
SaF_ST17 <- function(Sa, PGA, Tn, vs30, vref = 760) {
  
  stopifnot(length(Sa) == length(PGA))
  
  ## ---- special case: no amplification -----------------------------------
  if (vs30 == vref) {
    return(data.table(
      muLnSaF = log(Sa),   # SaF = Sa
      sdLnSaF = 0,
      ID      = "ST17"
    ))
  }
  
  ## ---- F(Tn,PGA,Vs30) statistics (original formulae) --------------------
  FModel <- F_ST17(PGA = PGA, Tn = Tn, vs30 = vs30, vref = vref)
  
  stopifnot(all(c("muLnF", "sdLnF") %in% names(FModel)),
            !any(is.na(FModel$sdLnF)))   # sd should be valid here
  
  ## ---- combine with Sa --------------------------------------------------
  muLnSaF <- log(Sa) + FModel$muLnF
  sdLnSaF <- FModel$sdLnF
  
  data.table(
    muLnSaF = muLnSaF,
    sdLnSaF = sdLnSaF,
    ID      = "ST17"
  )
}

# ============================================================
#  Stewart & Seyhan (2017) site-factor model  —  F_ST17()
# ============================================================
# nolint start
#' Non-linear site amplification factor **F<sub>ST17</sub>**
#'
#' Implements the period-dependent, Vs30-dependent model of Stewart & Seyhan (2017)
#' for the natural-log mean and standard deviation of the amplification factor
#' \eqn{F(T_n,\mathrm{PGA},V_{s30})}.
#'
#' *If \code{vs30 == vref}, the site factor equals 1, so the function returns
#'   \eqn{\mu_{\ln F}=0} and \code{sdLnF = NA} (to signal deterministic behaviour).*
#'
#' @param PGA   numeric **vector** – peak ground acceleration (in g).
#' @param Tn    numeric **scalar** – oscillator period (s).
#' @param vs30  numeric **scalar** – time-averaged shear-wave velocity to 30 m (m/s).
#' @param vref  numeric **scalar** – reference velocity (760 m/s or 3000 m/s), default 760.
#' @param Vl    numeric scalar – lower Vs30 bound for NL sigma interpolation (default 200).
#' @param Vu    numeric scalar – upper Vs30 bound for NL sigma interpolation (default 2000).
#'
#' @return A \link[data.table]{data.table} with columns:
#'   \describe{
#'     \item{\code{muLnF}}{natural-log mean of \eqn{F}}
#'     \item{\code{sdLnF}}{natural-log standard deviation of \eqn{F}
#'                         (set to \code{NA} when \code{vs30 == vref})}
#'     \item{\code{ID}}{character string \code{"ST17"}}}
#'
#' @section References:
#' Stewart, J.P. & Seyhan, E. (2017) “Semi-empirical nonlinear site
#' amplification model for global application.” *Earthquake Spectra* **33**(1).
#'
#' @export
F_ST17 <- function(PGA, Tn, vs30,
                   vref = 760,
                   Vl   = 200,
                   Vu   = 2000) {
  
  ## ----- basic input checks ---------------------------------------------
  ok <- length(vref) == 1L && length(vs30) == 1L && vref %in% c(760, 3000)
  stopifnot(ok)
  
  ## ----- deterministic case (vs30 == vref) -------------------------------
  if (vs30 == vref) {
    return(data.table(
      muLnF = 0,
      sdLnF = NA,
      ID    = "ST17"
    ))
  }
  
  ## ---------- interpolated coefficient tables (unchanged) ---------------
  cI     <- stats::approxfun(
    x = c(0.00, 0.08, 0.10, 0.20, 0.30, 0.40, 0.50, 0.80,
          1.00, 2.00, 3.00, 4.00, 5.00, 20.00),
    y = c(-0.359, -0.359, -0.344, -0.408, -0.416, -0.418, -0.446, -0.527,
          -0.554, -0.580, -0.593, -0.597, -0.582, -0.582),
    ties = mean
  )
  V1I    <- stats::approxfun(
    x = c(0.00, 0.08, 0.10, 0.20, 0.30, 0.40, 0.50, 0.80,
          1.00, 2.00, 3.00, 4.00, 5.00, 20.00),
    y = c(333, 333, 319, 314, 200, 200, 200, 230,
          278, 286, 313, 322, 325, 325),
    ties = mean
  )
  V2I    <- stats::approxfun(
    x = c(0.00, 0.08, 0.10, 0.20, 0.30, 0.40, 0.50, 0.80,
          1.00, 2.00, 3.00, 4.00, 5.00, 20.00),
    y = c(760, 760, 760, 760, 760, 760, 764, 942,
          1103, 1201, 876, 881, 855, 855),
    ties = mean
  )
  VfI    <- stats::approxfun(
    x = c(0.00, 0.08, 0.10, 0.20, 0.30, 0.40, 0.50, 0.80,
          1.00, 2.00, 3.00, 4.00, 5.00, 20.00),
    y = c(333, 333, 319, 314, 250, 250, 280, 280,
          300, 300, 313, 322, 325, 325),
    ties = mean
  )
  sdVcI  <- stats::approxfun(
    x = c(0.00, 0.08, 0.10, 0.20, 0.30, 0.40, 0.50, 0.80,
          1.00, 2.00, 3.00, 4.00, 5.00, 20.00),
    y = c(0.268, 0.268, 0.270, 0.251, 0.225, 0.225, 0.225, 0.225,
          0.225, 0.259, 0.306, 0.340, 0.340, 0.340),
    ties = mean
  )
  sdLI   <- stats::approxfun(
    x = c(0.00, 0.08, 0.10, 0.20, 0.30, 0.40, 0.50, 0.80,
          1.00, 2.00, 3.00, 4.00, 5.00, 20.00),
    y = c(0.281, 0.281, 0.263, 0.306, 0.276, 0.275, 0.311, 0.334,
          0.377, 0.548, 0.538, 0.435, 0.400, 0.400),
    ties = mean
  )
  sdUI   <- stats::approxfun(
    x = c(0.00, 0.08, 0.10, 0.20, 0.30, 0.40, 0.50, 0.80,
          1.00, 2.00, 3.00, 4.00, 5.00, 20.00),
    y = c(0.472, 0.472, 0.470, 0.334, 0.381, 0.381, 0.323, 0.308,
          0.361, 0.388, 0.551, 0.585, 0.587, 0.587),
    ties = mean
  )
  F760I  <- stats::approxfun(
    x = c(0.00, 0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017,
          0.019, 0.02, 0.022, 0.023, 0.025, 0.027, 0.028, 0.031,
          0.033, 0.035, 0.038, 0.04, 0.043, 0.046, 0.05, 0.053,
          0.057, 0.061, 0.066, 0.071, 0.076, 0.081, 0.087, 0.093,
          0.1, 0.107, 0.115, 0.123, 0.132, 0.142, 0.152, 0.163,
          0.175, 0.187, 0.201, 0.215, 0.231, 0.248, 0.266, 0.285,
          0.305, 0.327, 0.351, 0.376, 0.404, 0.433, 0.464, 0.498,
          0.534, 0.572, 0.614, 0.658, 0.705, 0.756, 0.811, 0.87,
          0.933, 1, 1.072, 1.15, 1.233, 1.322, 1.417, 1.52, 1.63,
          1.748, 1.874, 2.009, 2.154, 2.31, 2.477, 2.656, 2.848,
          3.054, 3.275, 3.511, 3.765, 4.037, 4.329, 4.642, 4.977,
          5.337, 5.722, 6.136, 6.579, 7.055, 7.565, 8.111, 8.697,
          9.326, 10, 20),
    y = c(0.185, 0.185, 0.185, 0.185, 0.185, 0.185, 0.185, 0.185,
          0.185, 0.185, 0.185, 0.185, 0.189, 0.195, 0.203, 0.212,
          0.224, 0.238, 0.252, 0.267, 0.283, 0.3, 0.318, 0.337,
          0.356, 0.377, 0.4, 0.425, 0.454, 0.475, 0.512, 0.558,
          0.613, 0.674, 0.73, 0.76, 0.759, 0.714, 0.647, 0.586,
          0.534, 0.488, 0.449, 0.419, 0.39, 0.362, 0.332, 0.301,
          0.278, 0.27, 0.262, 0.242, 0.224, 0.209, 0.197, 0.186,
          0.175, 0.166, 0.157, 0.15, 0.142, 0.135, 0.127, 0.12,
          0.111, 0.103, 0.095, 0.088, 0.083, 0.08, 0.078, 0.078,
          0.077, 0.077, 0.078, 0.079, 0.079, 0.078, 0.076, 0.075,
          0.074, 0.073, 0.073, 0.072, 0.07, 0.068, 0.066, 0.065,
          0.065, 0.064, 0.063, 0.061, 0.06, 0.058, 0.057, 0.056,
          0.056, 0.055, 0.055, 0.053, 0.053),
    ties = mean
  )
  sd760I <- stats::approxfun(
    x = c(0.00, 0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017,
          0.019, 0.02, 0.022, 0.023, 0.025, 0.027, 0.028, 0.031, 0.033,
          0.035, 0.038, 0.04, 0.043, 0.046, 0.05, 0.053, 0.057, 0.061,
          0.066, 0.071, 0.076, 0.081, 0.087, 0.093, 0.1, 0.107, 0.115,
          0.123, 0.132, 0.142, 0.152, 0.163, 0.175, 0.187, 0.201, 0.215,
          0.231, 0.248, 0.266, 0.285, 0.305, 0.327, 0.351, 0.376, 0.404,
          0.433, 0.464, 0.498, 0.534, 0.572, 0.614, 0.658, 0.705, 0.756,
          0.811, 0.87, 0.933, 1, 1.072, 1.15, 1.233, 1.322, 1.417, 1.52,
          1.63, 1.748, 1.874, 2.009, 2.154, 2.31, 2.477, 2.656, 2.848,
          3.054, 3.275, 3.511, 3.765, 4.037, 4.329, 4.642, 4.977, 5.337,
          5.722, 6.136, 6.579, 7.055, 7.565, 8.111, 8.697, 9.326, 10, 20),
    y = c(0.434, 0.434, 0.434, 0.434, 0.434, 0.434, 0.434, 0.434, 0.434,
          0.434, 0.434, 0.434, 0.432, 0.429, 0.422, 0.414, 0.404, 0.393,
          0.387, 0.387, 0.39, 0.39, 0.381, 0.363, 0.34, 0.32, 0.308, 0.306,
          0.312, 0.322, 0.335, 0.346, 0.357, 0.366, 0.372, 0.37, 0.351,
          0.323, 0.284, 0.253, 0.234, 0.222, 0.214, 0.214, 0.207, 0.195,
          0.177, 0.156, 0.141, 0.131, 0.124, 0.117, 0.115, 0.112, 0.113,
          0.111, 0.105, 0.104, 0.113, 0.123, 0.132, 0.139, 0.138, 0.133,
          0.13, 0.128, 0.124, 0.118, 0.112, 0.108, 0.11, 0.114, 0.12,
          0.122, 0.124, 0.124, 0.118, 0.113, 0.112, 0.111, 0.109, 0.108,
          0.111, 0.116, 0.12, 0.122, 0.12, 0.116, 0.112, 0.108, 0.104,
          0.1, 0.096, 0.091, 0.087, 0.082, 0.078, 0.074, 0.07, 0.069, 0.069),
    ties = mean
  )
  f3I    <- stats::approxfun(
    x = c(0.0, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1., 2., 3., 4., 5., 10., 20.),
    y = c(0.16249, 0.16249, 0.15083, 0.12815, 0.1307, 0.09414, 0.09888, 0.07357,
          0.04367, 0.00164, 0.00746, 0.00269, 0.00242, 0.05329, 0.05329),
    ties = mean
  )
  f4I    <- stats::approxfun(
    x = c(0.0, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1., 2., 3., 4., 5., 10., 20.),
    y = c(-0.50667, -0.50667, -0.44661, -0.30481, -0.22825, -0.11591, -0.07793,
          -0.01592, -0.00478, -0.00236, -0.00626, -0.00331, -0.00256, -0.00631,
          -0.00631),
    ties = mean
  )
  f5I    <- stats::approxfun(
    x = c(0.0, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1., 2., 3., 4., 5., 10., 20.),
    y = c(-0.00273, -0.00273, -0.00335, -0.00488, -0.00655, -0.00872, -0.01028,
          -0.01515, -0.01823, -0.01296, -0.01043, -0.01215, -0.01325,
          -0.01403, -0.01403),
    ties = mean
  )
  VcI    <- stats::approxfun(
    x = c(0.0, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1., 2., 3., 4., 5., 10., 20.),
    y = c(2990, 2990, 2990, 1533, 1152, 1018, 938, 832, 951, 879, 894, 875,
          856, 837, 837),
    ties = mean
  )
  sdcI   <- stats::approxfun(
    x = c(0.0, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1., 2., 3., 4., 5., 10., 20.),
    y = c(0.12, 0.12, 0.12, 0.12, 0.15, 0.15, 0.15, 0.1, 0.06, 0.04, 0.04, 0.03,
          0.02, 0.02, 0.02),
    ties = mean
  )
  
  ## ---------- mean term --------------------------------------------------
  muLnPGA <- log(PGA)
  C7603000 <- if (vref == 760) 2.275 else 1.0
  
  if (vs30 <= V1I(Tn)) {
    muL <- cI(Tn) * log(V1I(Tn) / 760)
  } else if (vs30 <= V2I(Tn)) {
    muL <- cI(Tn) * log(vs30 / 760)
  } else {
    muL <- cI(Tn) * log(V2I(Tn) / 760) + cI(Tn) / 2 * log(vs30 / V2I(Tn))
  }
  
  a <- f4I(Tn) * (exp(f5I(Tn) * (min(vs30, 3000) - 360)) -
                    exp(f5I(Tn) * (3000 - 360)))
  b <- C7603000 * f3I(Tn)
  
  muNL <- if (vs30 < VcI(Tn)) a * log(exp(muLnPGA) / b + 1) else 0
  muI  <- if (vref == 3000) F760I(Tn) else 0
  
  muLnF <- muL + muI + muNL
  
  ## ---------- sigma term -------------------------------------------------
  if (vs30 <= VfI(Tn)) {
    sdL <- sdLI(Tn) - 2 * (sdLI(Tn) - sdVcI(Tn)) *
      (vs30 - Vl) / (VfI(Tn) - Vl) +
      (sdLI(Tn) - sdVcI(Tn)) * (vs30 - Vl)^2 / (VfI(Tn) - Vl)^2
  } else if (vs30 <= V2I(Tn)) {
    sdL <- sdVcI(Tn)
  } else {
    sdL <- sdVcI(Tn) + (sdUI(Tn) - sdVcI(Tn)) *
      (vs30 - V2I(Tn))^2 / (Vu - V2I(Tn))^2
  }
  
  if (vs30 < 300) {
    sdNL <- sdcI(Tn)
  } else if (vs30 < 1000) {
    sdNL <- -sdcI(Tn) * log(vs30 / 300) / log(1000 / 300) + sdcI(Tn)
  } else {
    sdNL <- 0
  }
  
  sdI <- if (vref == 3000) sd760I(Tn) else 0
  sdLnF <- sqrt(sdL^2 + sdI^2 + sdNL^2)
  
  ## ---------- output -----------------------------------------------------
  data.table(
    muLnF = muLnF,
    sdLnF = sdLnF,
    ID    = "ST17"
  )
}

# nolint end


