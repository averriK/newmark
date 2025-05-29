# nolint start
Dn_model <- function(PGA,Sa,Tn, ky = NULL, Ts = NULL,  Mw = NULL, PGV = NULL, AI = NULL, kymin = 0.005, kymax = 0.5, n = 30) {
  if (is.null(ky)) {
    ky <- seq(from = log(kymin), to = log(kymax), length.out = n) |>
      exp() |>
      round(digits = 4) |>
      unique()
  }

  DT <- data.table()

  
  # DECOUPLED METHODS
  DT <- list(DT,Dn_AM88(PGA = PGA, ky = ky)) |> rbindlist(use.names = TRUE, fill = TRUE)
  DT <- list(DT,Dn_JB07(PGA = PGA, ky = ky,AI=AI)) |> rbindlist(use.names = TRUE, fill = TRUE)
  DT <- list(DT, Dn_SR08(PGA = PGA, ky = ky,AI=AI)) |> rbindlist(use.names = TRUE, fill = TRUE)
  DT <- list(DT,Dn_YG91(PGA = PGA, ky = ky)) |> rbindlist(use.names = TRUE, fill = TRUE)
  
  # COUPLED METHODS
  if(!is.null(Ts)){
    uhs <- data.table(Sa,Tn)
    #
    uhs <- list(uhs[Tn > 0] ,data.table(Tn = 0.01, Sa = PGA),data.table(Tn = 0.005, Sa = PGA)) |> rbindlist(use.names = TRUE, fill = TRUE)
    uhs <- uhs[order(Tn)]
    
    #
    shift <- 1.5
    Sa <- stats::approx(x = log(uhs$Tn), y = log(uhs$Sa), xout = log(shift*Ts))$y |> exp()
    
    DT <- list(DT,Dn_BT07(Ts = Ts, Sa = Sa, Mw = Mw, ky = ky)) |> rbindlist(use.names = TRUE, fill = TRUE)
    #
    shift <- 1.5
    Sa <- stats::approx(x = log(uhs$Tn), y = log(uhs$Sa), xout = log(shift*Ts))$y |> exp()
    
    DT <- list(DT,Dn_BM17(Ts = Ts, Sa = Sa, Mw = Mw, ky = ky)) |> rbindlist(use.names = TRUE, fill = TRUE) # Subduction)
    
    # 
    shift <- 1.3
    Sa <- stats::approx(x = log(uhs$Tn), y = log(uhs$Sa), xout = log(shift*Ts))$y |> exp()
    
    DT <- list(DT,Dn_BM19(Ts = Ts, Sa = Sa, PGA=PGA,Mw = Mw, ky = ky,PGV=PGV)) |> rbindlist(use.names = TRUE, fill = TRUE) # Shallow Crustal)
  }

  return(DT)
}

# ---------------------------------------------------------------------------
#  Bray & Macedo (2017) empirical displacement model – subduction events
#
#  Arguments
#    ky : yield acceleration ratio (dimensionless)
#    Sa : spectral acceleration at shifted period 1.5*Ts  [g]
#    Ts : fundamental period of slide mass  [s]
#    Mw : moment magnitude
#
#  Returns  data.table(muLnD, sdLnD, ID)
# ---------------------------------------------------------------------------
Dn_BM17 <- function(ky, Sa, Ts, Mw=7.5) {
  lnSa <- log(Sa)
  
  
  # --- period-dependent constants ------------------------------------------
  if (Ts < 0.1) {
    c1 <- -5.864
    c2 <- -9.421
    c3 <-   0.0      # short-period branch
  } else {
    c1 <- -6.896
    c2 <-  3.081
    c3 <-  -0.803    # long-period branch
  }
  # --- regression coefficients ---------------------------------------------
  a0 <- c1 + 0.550 * Mw + c2 * Ts + c3 * Ts^2 - 3.353 * log(ky) - 0.390 * log(ky)^2
  a1 <- 3.060 + 0.538 * log(ky)
  a2 <- -0.225
  # --- median & dispersion --------------------------------------------------
  muLnD <- a0 + a1 * lnSa + a2 * lnSa^2
  sdLnD <- 0.73
  data.table(muLnD = muLnD, sdLnD = sdLnD, ID = "BM17")
}

# ---------------------------------------------------------------------------
#  Bray & Macedo (2019) displacement model – shallow-crustal events
#
#  Arguments
#    ky  : yield-acceleration ratio
#    Ts  : fundamental period [s]
#    Sa  : spectral acceleration at shifted period 1.3 Ts [g]
#    Mw  : moment magnitude
#    PGA : peak ground acceleration [g]  (optional; used to back-calculate PGV)
#    PGV : peak ground velocity [cm/s] (optional; overrides PGA-based estimate)
#
#  Returns  data.table(muLnD, sdLnD, ID)
# ---------------------------------------------------------------------------
Dn_BM19 <- function(ky, Ts, Sa, PGA, Mw = 6.5, PGV) {
  # Validate vector lengths
  n <- length(Sa)
  if (!(length(PGA) == n && length(PGV) == n)) {
    stop("Sa, PGA, and PGV must have the same length.")
  }
  
  lnSa <- log(Sa)
  I <- (PGV <= 115) # Logical vector, TRUE=ordinary, FALSE=pulse
  
  # Coefficients (scalar, since Ts is scalar)
  c1 <- ifelse(Ts < 0.1, -4.551, -5.894)
  c2 <- ifelse(Ts < 0.1, -9.690, 3.152)
  c3 <- ifelse(Ts < 0.1, 0.000, -0.910)
  
  lnPGV_term <- ifelse(I, 0, 1.0)
  PGV_shift  <- ifelse(I, 0, -log(115))
  sdLnD      <- 0.74
  
  a0 <- c1 + 0.607 * Mw + c2 * Ts + c3 * Ts^2 -
    2.491 * log(ky) - 0.245 * log(ky)^2 +
    lnPGV_term * log(PGV) + PGV_shift
  
  a1 <- 2.703 + 0.344 * log(ky)
  a2 <- -0.089
  
  muLnD <- a0 + a1 * lnSa + a2 * lnSa^2
  
  data.table(
    muLnD = muLnD,
    sdLnD = sdLnD,
    ID = "BM19"
  )
}






# ---------------------------------------------------------------------------
#  Bray & Travasarou (2007) flexible sliding-block model
#
#  Arguments
#    ky : yield-acceleration ratio
#    Sa : spectral acceleration at shifted period 1.5 Ts [g]
#    Ts : fundamental period [s]    (used only to obtain Sa; not in formula)
#    Mw : moment magnitude
#
#  Returns  data.table(muLnD, sdLnD, ID)
# ---------------------------------------------------------------------------
Dn_BT07 <- function(ky, Sa,Ts, Mw=6.5) {
  
  lnSa <- log(Sa)
  lnky  <- log(ky)
  muLnD <- -1.10 - 2.83 * lnky - 0.333 * lnky^2 + 0.566 * lnky * lnSa + 3.04 * lnSa - 0.244 * lnSa^2 +  0.278 * (Mw - 7)
  sdLnD <- 0.66
  data.table(muLnD = muLnD, sdLnD = sdLnD, ID = "BT07")
}

# ---------------------------------------------------------------------------
#  Jibson (2007) empirical displacement model
#
#  Arguments
#    PGA : peak ground acceleration [g]
#    AI  : Arias intensity [m/s] (optional; computed from PGA if missing)
#    ky  : yield acceleration [g]
#
#  Returns  data.table(muLnD, sdLnD, ID)
# ---------------------------------------------------------------------------
Dn_JB07 <- function(PGA, ky, AI = NULL) {
  if (is.null(AI)) {
    AI <- (PGA^1.9228) * exp(2.6109)           # AI back-calculated [m/s]
  }
  r        <- ky / PGA
  muLog10D <- 0.561 * log10(AI) - 3.833 * log10(r) - 1.474
  sdLog10D <- 0.616
  data.table(
    muLnD = muLog10D * log(10),
    sdLnD = sdLog10D * log(10),
    ID    = "JB07"
  )
}

# ---------------------------------------------------------------------------
#  Saygili & Rathje (2008) displacement model
#
#  Arguments
#    PGA : peak ground acceleration [g]
#    ky  : yield acceleration [g]
#    AI  : Arias intensity [m/s] (optional; computed from PGA if missing)
#
#  Returns  data.table(muLnD, sdLnD, ID)
# ---------------------------------------------------------------------------
Dn_SR08 <- function(PGA, ky, AI = NULL) {
  if (is.null(AI)) {
    AI <- (PGA^1.9228) * exp(2.6109)           # AI back-calculated [m/s]
  }
  r      <- ky / PGA
  muLnD  <- 2.39 - 5.24 * r - 18.78 * r^2 + 42.01 * r^3 - 29.15 * r^4 -
    1.56 * log(PGA) + 1.38 * log(AI)
  sdLnD  <- 0.46 + 0.56 * r
  data.table(muLnD = muLnD, sdLnD = sdLnD, ID = "SR08")
}

# ---------------------------------------------------------------------------
#  Yegian et al. (1991) rigid sliding-block model
#
#  Arguments
#    PGA : peak ground acceleration [g]
#    ky  : yield acceleration [g]
#
#  Returns  data.table(muLnD, sdLnD, ID)
# ---------------------------------------------------------------------------
Dn_YG91 <- function(PGA, ky) {
  r <- ky / PGA
  muLog10D <- 0.22 - 10.12 * r + 16.38 * r^2 - 11.48 * r^3
  sdLog10D <- 0.45
  data.table(
    muLnD = muLog10D * log(10),
    sdLnD = sdLog10D * log(10),
    ID    = "YG91"
  )
}


# ---------------------------------------------------------------------------
#  Ambraseys & Menu (1988) rigid sliding-block displacement model
#
#  Arguments
#    PGA : peak ground acceleration  [g]
#    ky  : yield acceleration        [g]   ( = k_y * g )
#
#  Returns
#    data.table with
#       muLnD : mean of ln D   [–]   (D in centimetres)
#       sdLnD : σ_lnD         [–]
#       ID    : model tag     [chr]
# ---------------------------------------------------------------------------
Dn_AM88 <- function(PGA, ky) {
  r  <- ky / PGA                       # dimensionless strength ratio
  # --- model coefficients (Table 2, Ambraseys & Menu 1988) ------------------
  a1 <- 0.90
  a2 <- 2.53
  a3 <- -1.09
  # --- median in log10 space -----------------------------------------------
  muLog10D <- a1 + log10( (1 - r)^a2 * r^a3 )
  sdLog10D <- 0.30                    # σ_log10 D
  # --- convert to natural-log ----------------------------------------------
  muLnD <- muLog10D * log(10)
  sdLnD <- sdLog10D * log(10)
  data.table(muLnD = muLnD, sdLnD = sdLnD, ID = "AM88")
}


