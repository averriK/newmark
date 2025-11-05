#' Baker and Jayaram (2009) intra-record spectral correlation model
#'
#' Computes the Pearson correlation coefficient between ln(Sa) at two periods
#' (typically between PGA and Sa(Tn)) using the exponential decay model proposed
#' by Jayaram and Baker (2009), calibrated as a function of rupture distance.
#'
#' The formula used is:
#'   rho(T1, T2) = exp(-alpha * abs(log(T1 / T2)))
#' where alpha depends on Rrup as follows:
#'   - Rrup < 10 km     -> alpha = 0.35
#'   - 10 <= Rrup <= 20 -> alpha = 0.30
#'   - Rrup > 20 km     -> alpha = 0.25
#'
#' Reference: Jayaram, N., & Baker, J. W. (2009). Correlation model for spatially
#' distributed ground-motion intensities. Earthquake Engineering & Structural
#' Dynamics, 38(8), 951–972. https://doi.org/10.1002/eqe.876
#'
#' @param Tn Numeric scalar. Target spectral period in seconds.
#' @param T0 Numeric scalar. Reference period. Default is 0.01 (PGA).
#' @param Rrup Numeric scalar. Rupture distance (km). Determines alpha.
#'
#' @return Numeric scalar correlation coefficient rho(T0, Tn)
#'
#' @export
rhoBJ <- function(Tn, T0 = 0.01, Rrup = 100) {
  stopifnot(length(Tn) == 1L, length(T0) == 1L, length(Rrup) == 1L)
  alpha <- if (Rrup < 10) {
    0.35
  } else if (Rrup <= 20) {
    0.30
  } else {
    0.25
  }

  exp(-alpha * abs(log(T0 / Tn)))
}
