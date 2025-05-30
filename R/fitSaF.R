# ============================================================
#  Public wrapper: fitSaF()
#  (mirrors fitDn() workflow)
# ============================================================
# nolint start
#' Quantiles of SaF = Sa × F(ST17) at period Tn
#'
#' @param UHS  data.table – uniform-hazard spectrum (Tn, Sa, p).
#' @param vs30 numeric scalar – site Vs30 (m/s).
#' @param Tn   numeric scalar – oscillator period (s).
#' @param NS   integer ≥ 1 – Monte-Carlo samples (default 30).
#' @param vref numeric scalar – reference Vs (default 760 m/s).
#'
#' @return data.table(p, SaF) with probability labels plus row p == "mean".
#' @export
fitSaF <- function(UHS,
                   vs30,
                   Tn,
                   NS   = 30,
                   vref = 760) {
  
  ## ---------- 1. sample Sa(Tn) and PGA ----------------------------------
  Sa_Tn <- sampleSa(UHS, Td = Tn,   n = NS)$Sa
  PGA   <- sampleSa(UHS, Td = 0.01, n = NS)$Sa  # 0.01 s ≈ PGA
  
  ## ---------- 2. SaF Monte-Carlo sampling -------------------------------
  SaFTable <- getSaF(
    SaF_ST17,
    Sa   = Sa_Tn,
    PGA  = PGA,
    Tn   = Tn,
    vs30 = vs30,
    vref = vref,
    n    = NS
  )
  
  ## ---------- 3. empirical quantiles ------------------------------------
  probs <- UHS[p != "mean", unique(as.numeric(p))]
  
  SaF_q <- data.table(
    p   = probs,
    SaF = stats::quantile(
      SaFTable$SaF,
      probs = probs,
      names = FALSE,
      type  = 7)
  )
  
  SaF_mean <- data.table(p = "mean",
                         SaF = mean(SaFTable$SaF))
  
  rbindlist(list(SaF_q, SaF_mean))
}

# ============================================================
#  Internal helper: getSaF()
#  (signature & style identical to getDn())
# ============================================================
# nolint start
#' Draw Monte-Carlo samples of SaF
#' @keywords internal
getSaF <- function(.fun, ..., n = 1) {
  
  SaFModel <- .fun(...)                     # identical idiom to getDn()
  stopifnot(all(c("muLnSaF", "sdLnSaF", "ID") %in% names(SaFModel)))
  
  SaFModel[
    , .(ID, LnSaF = rnorm(n, muLnSaF, sdLnSaF)),
    by = .I
  ][
    , .(sample = .I, ID, SaF = exp(LnSaF))]
}



# nolint end
