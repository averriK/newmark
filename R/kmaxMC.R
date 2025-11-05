#' Probabilistic k_max at displacement targets (Monte Carlo)
#'
#' Builds a sample of Dn(ky) from the provided partial quantiles and inverts
#' Dn(ky)=Da to obtain the distribution of k_max for each target Da.
#' For each simulation, one Uniform(0,1) is drawn and reused across all ky.
#' For every ky, ln-quantiles are reconstructed with buildQSpline() from the
#' numeric p-rows, then rescaled to match the provided "mean" at that ky.
#' The simulated Dn(ky) is forced non-increasing in ky before inversion.
#'
#' @param Dn       numeric, displacements paired with ky and p (same length)
#' @param ky       numeric, yield acceleration (same length as Dn)
#' @param p        character/numeric labels; must include numeric quantiles and "mean"
#' @param Da       numeric vector of displacement targets (e.g. c(2.5,25,250))
#' @param p_TARGET numeric/character vector of probability levels to report
#'                 (e.g. c(0.16,"mean",0.84))
#' @param ns       integer, Monte Carlo size (default 5e4)
#'
#' @return data.table with columns: Da, p, kmax
#' @import data.table
#' @export
kmaxMC <- function(Dn, ky, p,
                   Da       = c(2.5, 25, 250),
                   p_TARGET = c(0.16, "mean", 0.84),
                   ns       = 5e4) {

  KY <- sort(unique(ky))
  if (length(KY) < 2L) {
    return(data.table::data.table(
      Da = rep(Da, each = length(p_TARGET)),
      p  = rep(as.character(p_TARGET), times = length(Da)),
      kmax = NA_real_
    ))
  }

  J   <- split(seq_along(ky), ky)   # row indices by ky
  Q   <- vector("list", length(KY)) # u -> ln(Dn) per ky
  mu  <- rep(NA_real_, length(KY))  # mean Dn per ky ("mean" row)

  for (j in seq_along(KY)) {
    idx <- J[[ as.character(KY[j]) ]]
    pj  <- p[idx]
    qj  <- Dn[idx]

    im <- which(pj == "mean" & is.finite(qj) & qj > 0)
    if (length(im) == 1L) mu[j] <- qj[im]

    pn <- suppressWarnings(as.numeric(pj))
    ok <- which(is.finite(pn) & pn > 0 & pn < 1 &
                is.finite(qj) & qj > 0)
    if (length(ok) < 2L) {
      return(data.table::data.table(
        Da = rep(Da, each = length(p_TARGET)),
        p  = rep(as.character(p_TARGET), times = length(Da)),
        kmax = NA_real_
      ))
    }
    DQ <- data.table::data.table(p = pn[ok], Sa = qj[ok])[order(p)]
    Q[[j]] <- buildQSpline(DQ)  # empirical inverse: u in (0,1) -> ln(Dn)
  }

  U <- stats::runif(ns)
  D <- matrix(NA_real_, nrow = ns, ncol = length(KY))
  for (j in seq_along(KY)) {
    x <- exp(Q[[j]](U))
    if (is.finite(mu[j])) {
      m <- mean(x)
      if (is.finite(m) && m > 0) x <- x * (mu[j] / m) # re-center to "mean"
    }
    D[, j] <- x
  }

  ord <- order(KY)
  kY  <- KY[ord]
  K   <- matrix(NA_real_, nrow = ns, ncol = length(Da))

  for (s in seq_len(ns)) {
    z <- D[s, ord]
    z <- cummin(z)          # enforce Dn non-increasing in ky
    xr <- rev(z)            # make x increasing for approx()
    yr <- rev(kY)
    for (k in seq_along(Da)) {
      K[s, k] <- stats::approx(xr, yr, xout = Da[k],
                               rule = 2L, ties = "ordered")$y
    }
  }

  pc <- as.character(p_TARGET)
  pn <- suppressWarnings(as.numeric(pc))
  inum <- which(is.finite(pn) & pn > 0 & pn < 1)
  hasm <- any(pc == "mean")

  OUT <- vector("list", length(Da))
  for (k in seq_along(Da)) {
    dt <- data.table::data.table(Da = Da[k], p = pc, kmax = NA_real_)
    if (length(inum)) {
      qv <- stats::quantile(K[, k], probs = pn[inum], na.rm = TRUE, names = FALSE)
      dt$kmax[inum] <- as.numeric(qv)
    }
    if (hasm) {
      i <- which(pc == "mean")
      dt$kmax[i] <- mean(K[, k], na.rm = TRUE)
    }
    OUT[[k]] <- dt
  }

  OUT <- data.table::rbindlist(OUT, use.names = TRUE)
  data.table::setorder(OUT, Da, p)
  OUT[]
}
