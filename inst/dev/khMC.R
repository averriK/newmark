#' Probabilistic k_h = k_max / PGA at displacement targets (Monte Carlo)
#'
#' Builds a sample of Dn(ky) from the provided partial quantiles and a sample
#' of PGA from its quantiles. For each Monte Carlo draw:
#'   1) one U~Unif(0,1) is reused across all ky to sample Dn(ky) coherently;
#'   2) an independent V~Unif(0,1) samples PGA;
#'   3) invert Dn(ky)=Da to obtain k_max, then compute k_h = k_max / PGA.
#'
#' For every ky, ln-quantiles are reconstructed with buildQSpline() from the
#' numeric p-rows, then rescaled to match the provided "mean" at that ky.
#' The simulated Dn(ky) is forced non-increasing in ky before inversion.
#' PGA is sampled from its own p–PGA grid (unique p rows), also re-centered
#' to its "mean" if present.
#'
#' @param Dn       numeric, displacements paired with ky and p (same length)
#' @param ky       numeric, yield acceleration (same length as Dn)
#' @param p        character/numeric labels aligned with Dn and PGA
#' @param PGA      numeric, PGA values aligned with p (often repeated across ky)
#' @param Da       numeric vector of displacement targets (e.g. c(2.5,25,250))
#' @param p_TARGET numeric/character vector of probability levels to report
#'                 (e.g. c(0.16,"mean",0.84))
#' @param ns       integer, Monte Carlo size (default 5e4)
#'
#' @return data.table with columns: Da, p, kh
#' @import data.table
#' @export
khMC <- function(Dn, ky, p, PGA,
                 Da       = c(2.5, 25, 250),
                 p_TARGET = c(0.16, "mean", 0.84),
                 ns       = 5e4) {

  ## ---- set up Dn(ky) splines (per ky) ---------------------------------
  KY <- sort(unique(ky))
  if (length(KY) < 2L) {
    return(data.table::data.table(
      Da = rep(Da, each = length(p_TARGET)),
      p  = rep(as.character(p_TARGET), times = length(Da)),
      kh = NA_real_
    ))
  }

  J   <- split(seq_along(ky), ky)
  Q   <- vector("list", length(KY))   # u -> ln(Dn) per ky
  mu  <- rep(NA_real_, length(KY))    # mean Dn per ky

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
        kh = NA_real_
      ))
    }
    DQ <- data.table::data.table(p = pn[ok], Sa = qj[ok])[order(p)]
    Q[[j]] <- buildQSpline(DQ)  # empirical inverse: u in (0,1) -> ln(Dn)
  }

  ## ---- sample Dn(ky) with one shared U across ky ----------------------
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

  ## ---- enforce Dn non-increasing in ky --------------------------------
  ord <- order(KY)
  kY  <- KY[ord]

  for (s in seq_len(ns)) {
    z <- D[s, ord]
    D[s, ord] <- cummin(z)
  }

  ## ---- PGA quantile spline (unique p rows) ----------------------------
  pn_all <- suppressWarnings(as.numeric(p))
  okp <- which(is.finite(pn_all) & pn_all > 0 & pn_all < 1 &
               is.finite(PGA)    & PGA    > 0)
  if (length(okp) < 2L) {
    return(data.table::data.table(
      Da = rep(Da, each = length(p_TARGET)),
      p  = rep(as.character(p_TARGET), times = length(Da)),
      kh = NA_real_
    ))
  }
  Ptab <- data.table::data.table(p = pn_all[okp], Sa = PGA[okp])[order(p)]
  Ptab <- Ptab[, .SD[1L], by = p]  # one row per p
  Qpga <- buildQSpline(Ptab)

  muP <- NA_real_
  imP <- which(p == "mean" & is.finite(PGA) & PGA > 0)
  if (length(imP)) muP <- unique(PGA[imP])[1L]

  ## ---- sample PGA independently ---------------------------------------
  V <- stats::runif(ns)
  Pg <- exp(Qpga(V))
  if (is.finite(muP)) {
    mP <- mean(Pg)
    if (is.finite(mP) && mP > 0) Pg <- Pg * (muP / mP)
  }

  ## ---- invert to k_max, then divide by PGA ----------------------------
  K  <- matrix(NA_real_, nrow = ns, ncol = length(Da))
  Kh <- matrix(NA_real_, nrow = ns, ncol = length(Da))

  for (s in seq_len(ns)) {
    z  <- D[s, ord]
    xr <- rev(z)          # increasing in Dn
    yr <- rev(kY)
    for (k in seq_along(Da)) {
      km <- stats::approx(xr, yr, xout = Da[k],
                          rule = 2L, ties = "ordered")$y
      K[s, k]  <- km
      Kh[s, k] <- km / Pg[s]
    }
  }

  ## ---- quantiles for kh -----------------------------------------------
  pc <- as.character(p_TARGET)
  pn <- suppressWarnings(as.numeric(pc))
  inum <- which(is.finite(pn) & pn > 0 & pn < 1)
  hasm <- any(pc == "mean")

  OUT <- vector("list", length(Da))
  for (k in seq_along(Da)) {
    dt <- data.table::data.table(Da = Da[k], p = pc, kh = NA_real_)
    if (length(inum)) {
      qv <- stats::quantile(Kh[, k], probs = pn[inum], na.rm = TRUE, names = FALSE)
      dt$kh[inum] <- as.numeric(qv)
    }
    if (hasm) {
      i <- which(pc == "mean")
      dt$kh[i] <- mean(Kh[, k], na.rm = TRUE)
    }
    OUT[[k]] <- dt
  }

  OUT <- data.table::rbindlist(OUT, use.names = TRUE)
  data.table::setorder(OUT, Da, p)
  OUT[]
}
