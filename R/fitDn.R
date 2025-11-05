###############################################################################
# fitDn.R ― Displacement model ensemble with Monte Carlo & epistemic weights
###############################################################################
#' Probabilistic Newmark displacements (Monte Carlo + epistemic mixing)
#'
#' @param uhs  data.table ― columns Tn, p, Sa  (must include p == "mean")
#' @param ky   numeric scalar, yield acceleration (g)
#' @param Ts   numeric scalar, fundamental period of sliding mass (s)
#' @param Mw   numeric scalar, moment magnitude (default 6.5)
#' @param NS   integer, number of Monte‑Carlo samples
#' @param Rrup numeric scalar, rupture distance for rhoBJ()
#'
#' @importFrom stats rnorm
#' @return data.table with columns p, Dn  (displacement in cm)
#' @export
fitDn <- function(uhs,
                  ky,
                  Ts,
                  Mw   = 6.5,
                  NS   = 30L,
                  Rrup = 100) {

  # ---------------------------------------------------------------------
  # 1. Sanear la entrada: sólo columnas canónicas
  # ---------------------------------------------------------------------
  uhs <- uhs[, .(Tn, p, Sa)]

  # ---------------------------------------------------------------------
  # 2. Periodos objetivo (PGA + tres Sd desplazadas)
  # ---------------------------------------------------------------------
  Td0  <- 0.00
  Td10 <- 1.0 * Ts
  Td13 <- 1.3 * Ts
  Td15 <- 1.5 * Ts
  Td   <- c(Td0, Td10, Td13, Td15)

  # ---------------------------------------------------------------------
  # 3. Interpolar (log‑log) Sa(p) donde falte
  # ---------------------------------------------------------------------
  SaTable <- rbindlist(lapply(Td, function(Ti) {
    if (Ti %in% uhs$Tn) {
      uhs[Tn == Ti]
    } else {
      interpolateSaTable(uhs, Tn = Ti)
    }
  }))

  # ---------------------------------------------------------------------
  # 4. Muestrear espectros correlacionados (Baker & Jayaram 2009)
  # ---------------------------------------------------------------------
  rho      <- sapply(Td[-1], rhoBJ, T0 = Td0, Rrup = Rrup)
  SaMatrix <- sampleSaCorr(SaTable, Tn = Td, rho = rho, NS = NS)

  PGA  <- SaMatrix[, 1L]
  Sa10 <- SaMatrix[, 2L]
  Sa13 <- SaMatrix[, 3L]
  Sa15 <- SaMatrix[, 4L]
  PGV  <- (PGA^1.0529) * exp(0.1241) * 100          # cm/s
  AI   <- (PGA^1.9228) * exp(2.6109)                # Arias intensity (m/s)

  # ---------------------------------------------------------------------
  # 5. Pesos epistemológicos de los modelos
  # ---------------------------------------------------------------------
  weight <- data.table(
    ID = c("YG91", "AM88", "JB07", "BT07", "SR08", "BM17", "BM19"),
    w  = c(1, 2, 2, 3, 3, 4, 4)
  )
  weight[, w := w / sum(w)]

  # ---------------------------------------------------------------------
  # 6. Runner por modelo con salvaguarda
  # ---------------------------------------------------------------------
  drawDn <- function(ID) {
    spec <- tryCatch({
      switch(ID,
             AM88 = Dn_AM88(PGA, ky),
             YG91 = Dn_YG91(PGA, ky),
             JB07 = Dn_JB07(PGA, ky, AI = AI),
             SR08 = Dn_SR08(PGA, ky, AI = AI),
             BT07 = Dn_BT07(ky = ky, Sa = Sa15, Ts = Ts, Mw = Mw),
             BM17 = Dn_BM17(ky = ky, Sa = Sa15, Ts = Ts, Mw = Mw),
             BM19 = Dn_BM19(ky = ky, Ts = Ts, Sa = Sa13,
                            PGA = PGA, PGV = PGV, Mw = Mw),
             stop("Unknown model ID: ", ID))
    }, error = function(e) NULL)

    # Aceptar vectores mu/sd de longitud NS
    if (!is.null(spec) &&
        all(is.finite(spec$muLnD)) &&
        all(is.finite(spec$sdLnD))) {

      data.table(ID = ID,
                 Dn = exp(rnorm(NS, spec$muLnD, spec$sdLnD)))
    } else {
      NULL
    }
  }

  # ---------------------------------------------------------------------
  # 7. Compilar salidas válidas
  # ---------------------------------------------------------------------
  DnList <- lapply(weight$ID, drawDn)
  DnList <- Filter(Negate(is.null), DnList)

  if (length(DnList) == 0L) {
    warning("fitDn(): all models failed for ky = ", ky)
    return(data.table(p = character(), Dn = numeric()))
  }

  DnTable <- rbindlist(DnList, use.names = TRUE)

  # --- ENLACE LIMPIO DE PESOS (evita reciclado) -------------------------
  DnTable <- merge(DnTable, weight, by = "ID", all.x = TRUE, sort = FALSE)

  # ---------------------------------------------------------------------
  # 8. Cuantiles y media ponderada
  # ---------------------------------------------------------------------
  p_vals <- sort(unique(as.numeric(uhs[p != "mean", p])))

  Dn_q <- data.table(
    p  = p_vals,
    Dn = Hmisc::wtd.quantile(DnTable$Dn,
                             weights = DnTable$w,
                             probs   = p_vals,
                             type    = "quantile",
                             na.rm   = TRUE)
  )

  Dn_mean <- data.table(
    p  = "mean",
    Dn = Hmisc::wtd.mean(DnTable$Dn,
                         weights = DnTable$w,
                         na.rm   = TRUE)
  )

  rbindlist(list(Dn_q, Dn_mean), use.names = TRUE)
}
