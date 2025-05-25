Dn_BM17 <- function(  ky,Sa,Tn, Mw) {
  # Bray Macedo Model 2017
  # SUBDUCTION EARTHQUAKE
  
  # Ts <- 1.5 * Tn 
  # Ts is not used in this model!
  # PGA in [g]
  # Dn: [cm]
  muLnSa <- log(Sa) # Sa(Ts)
  if (Tn < 0.1) {
    # Model 1
    c1 <- -5.864
    c2 <- -9.421
    c3 <- 0.0
  } else if (Tn >= 0.1) {
    # Model 2
    c1 <- -6.896
    c2 <- 3.081
    c3 <- -0.803
  }
  a0 <- c1 + 0.550 * Mw + c2 * Tn + c3 * Tn^2 - 3.353 * log(ky) - 0.390 * log(ky)^2
  a1 <- 3.060 + 0.538 * log(ky)
  a2 <- -0.225

  # MODEL ----
  muLnD <- a0 + a1 * muLnSa + a2 * muLnSa^2
  sdLnD <- 0.73
  return( data.table(muLnD=muLnD,sdLnD=sdLnD,ID="BM17"))
}
