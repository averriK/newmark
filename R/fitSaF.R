#' Site-amplified AUXtral acceleration (Stewart & Seyhan 2017)
#'
#' Computes site-amplified AUXtral acceleration using two modes:
#' (A) If the input UHS contains only p == "mean", rock Sa is treated as
#'     deterministic and the dispersion in SaF comes solely from F_ST17
#'     (sdLnF), sampled on the natural-log scale.
#' (B) If the input UHS contains numeric quantiles (e.g., "0.10", "0.50", ...),
#'     ln Sa(Tn) is sampled from those quantiles using an empirical inverse
#'     quantile function (buildQSpline). Dependence between PGA and Sa(Tn)
#'     is represented by a Gaussian "star" copula (Baker & Jayaram, 2009)
#'     with correlation rhoBJ(Tn, T0 = 0, Rrup). For each Monte Carlo draw,
#'     F_ST17 is evaluated with the paired PGA draw, and ln SaF = ln Sa + ln F
#'     with ln F ~ Normal(muLnF, sdLnF). Independence between ln Sa and ln F
#'     is assumed.
#'
#' Requirements:
#' - The input table must have columns Tn, p, Sa.
#' - It must include exactly one row (Tn == 0, p == "mean") providing PGA on rock.
#'
#' Output columns: Tn, p, Sa, SaF, AF, where AF = SaF / Sa.
#' In mode (B), the same numeric p values found in the input are returned plus
#' "mean". In mode (A), the set of p values is taken from p_TARGET plus "mean".
#'
#' @param uhs       data.table with columns Tn, p, Sa.
#' @param vs30      target Vs30 (m/s).
#' @param vref      reference Vs30 (m/s), default 760.
#' @param ns        Monte Carlo size (default 1000).
#' @param models    kept for signature compatibility (unused).
#' @param Rrup      rupture distance (km) for rhoBJ() in the quantile case (default 100).
#' @param p_TARGET  probabilities to report when only 'mean' is provided,
#'                  default c(0.05, 0.10, 0.16, 0.50, 0.84, 0.90, 0.95).
#'
#' @return data.table with columns Tn, p, Sa, SaF, AF.
#' @export
fitSaF <- function(uhs,
                       vs30,
                       vref = 760,
                       ns = 1000,
                       models = "ST17",
                       Rrup = 100,
                       p_TARGET = c(0.05,0.10, 0.16,0.50,0.84, 0.90, 0.95)) {

  stopifnot(data.table::is.data.table(uhs))
  if (!all(c("Tn","p","Sa") %in% names(uhs))) {
    stop("fitSaF: 'uhs' must have columns Tn, p, Sa")
  }

  ## PGA (roca, media única)
  PGA <- uhs[Tn == 0 & p == "mean", Sa]
  if (length(PGA) != 1L) {
    stop("fitSaF: need exactly one (Tn == 0, p == 'mean') row for PGA")
  }

  ## Detección del tipo de entrada (ramificación única)
  p_num_all <- suppressWarnings(as.numeric(uhs$p))
  have_quantiles <- any(!is.na(p_num_all) & p_num_all > 0 & p_num_all < 1)

  Tn <- sort(unique(uhs$Tn))
  OUT <- data.table::data.table(Tn = numeric(0), p = character(0),
                                Sa = numeric(0), SaF = numeric(0))

  if (!have_quantiles) {
    # ---------------------- (A) sólo p == "mean" -----------------------
    for (t in Tn) {
      Sa <- uhs[Tn == t & p == "mean", Sa]
      if (length(Sa) != 1L) { next } # Manejar Tns sin datos

      AUX  <- F_ST17(PGA = PGA, Tn = t, vs30 = vs30, vref = vref)
      muLnF <- AUX$muLnF
      sdLnF <- AUX$sdLnF

      X <- stats::rnorm(ns, mean = muLnF, sd = sdLnF)  # ln F
      Z <- log(Sa) + X                                 # ln SaF

      OUT <- data.table::rbindlist(list(
        OUT,
        data.table::data.table(
          Tn  = t,
          p   = sprintf("%.2f", p_TARGET),
          Sa  = Sa,
          SaF = as.numeric(stats::quantile(exp(Z), p_TARGET, type = 7))
        ),
        data.table::data.table(
          Tn  = t,
          p   = "mean",
          Sa  = Sa,
          SaF = mean(exp(Z))
        )
      ), use.names = TRUE)
    }

  } else {
    # ---------------------- (B) hay cuantiles numéricos -----------------
    ## Preparar el muestreo de PGA una única vez
    DP0 <- uhs[Tn == 0]
    p0  <- suppressWarnings(as.numeric(DP0$p))
    id0 <- which(!is.na(p0) & p0 > 0 & p0 < 1)
    DP0q <- DP0[id0, .(p = as.numeric(p), Sa)]
    DP0q <- DP0q[order(p)]
    LnQ0 <- buildQSpline(DP0q)  # u -> ln PGA

    for (t in Tn) {
      DT <- uhs[Tn == t]
      Sa <- DT[p == "mean", Sa]
      if (length(Sa) != 1L) { next } # Manejar Tns sin datos

      ## inverso empírico en ln‑espacio para Sa(Tn)
      p_num <- suppressWarnings(as.numeric(DT$p))
      idx   <- which(!is.na(p_num) & p_num > 0 & p_num < 1)
      DTn <- DT[idx, .(p = as.numeric(p), Sa)]
      DTn <- DTn[order(p)]
      LnQn <- buildQSpline(DTn)  # u in (0,1) -> ln Sa(Tn)

      ## cópula gaussiana “estrella” PGA–Sa(Tn)
      rho <- rhoBJ(Tn = t, T0 = 0, Rrup = Rrup)
      Sig <- matrix(c(1, rho, rho, 1), nrow = 2L)
      Z2  <- mvtnorm::rmvnorm(ns, sigma = Sig)  # (Z0, Zn)
      U0  <- stats::pnorm(Z2[, 1L])
      U1  <- stats::pnorm(Z2[, 2L])

      ## mapeo a ln‑marginales y recentering
      lnPGA_sample <- LnQ0(U0)
      lnPGA_sample <- log(exp(lnPGA_sample) * (PGA / mean(exp(lnPGA_sample))))

      lnSa_sample  <- LnQn(U1)
      lnSa_sample  <- log(exp(lnSa_sample) * (Sa / mean(exp(lnSa_sample))))

      ## ST17 con PGA por muestreo
      AUX_vec <- F_ST17(PGA = exp(lnPGA_sample), Tn = t, vs30 = vs30, vref = vref)
      lnF_sample   <- stats::rnorm(ns, mean = AUX_vec$muLnF, sd = AUX_vec$sdLnF)
      lnSaF <- lnSa_sample + lnF_sample

      ## mismos p que trae el hazard en este Tn + media
      OUT <- data.table::rbindlist(list(
        OUT,
        data.table::data.table(
          Tn  = t,
          p   = DT$p[idx],
          Sa  = DT$Sa[idx],
          SaF = as.numeric(stats::quantile(exp(lnSaF), p_num[idx], type = 7))
        ),
        data.table::data.table(
          Tn  = t,
          p   = "mean",
          Sa  = Sa,
          SaF = mean(exp(lnSaF))
        )
      ), use.names = TRUE)
    }
  }

  OUT[, AF := SaF / Sa]
  data.table::setorder(OUT, Tn, p)
  OUT[]
}
