Dn_BM19 <- function(ky,Tn, Sa,Mw,PGA=NULL,PGV=NULL) {
  # Bray Macedo Model 2019 ----
  # SHALOW CRUSTAL EARTHQUAKE
  # PGV in [cm/s]
  # PGA in [g]
  # Dn: [cm]
  # Ts <- 1.3 * Tn 
  # Ts is not used in this model
  muLnSa <- log(Sa) #Sa(Ts)
  
  if(is.null(PGV) & !is.null(PGA) & PGA>0){
    b0 <- 0.1241
    b1 <- 1.0529
    b2 <- 0.0
    PGV <- ((PGA^b1) * exp(b0 + (b2 / PGA)))*100 # [cm/s]
    # PGV <- 55 * PGA # Idriss [cm/s]
    
  }
  
  
  if (PGV <= 115) {
    # Ordinary GM
    if (Tn < 0.1) {
      # Model 1
      c1 <- -4.551
      c2 <- -9.690
      c3 <- 0.0
    } else {
      # Model 2
      c1 <- -5.894
      c2 <- 3.152
      c3 <- -0.910
    }
    c4 <- +0.0
    c5 <- +0.0
    sdLnDo <- 0.74
  }
  if (PGV > 115) {
    if (Tn < 0.1) {
      # Model 1
      c1 <- -4.551
      c2 <- -9.690
      c3 <- +0.0
    } else if (Tn >= 0.1) {
      # Model 2
      c1 <- -5.894
      c2 <- 3.152
      c3 <- -0.910
    } else {
      stop()
    }
    c4 <- +1.0
    c5 <- -4.75
    sdLnDo <- 0.74
  }
  
  a0 <- c1 + 0.607 * Mw + c2 * Tn + c3 * Tn^2 - 2.491 * log(ky) - 0.245 * log(ky)^2 + c4 * log(PGV) + c5
  a1 <- 2.703 + 0.344 * log(ky)
  a2 <- -0.089
  
  muLnD <- a0 + a1 * muLnSa + a2 * muLnSa^2
  sdLnD <- sdLnDo
  return( data.table(muLnD=muLnD,sdLnD=sdLnD,ID="BM19"))
}
