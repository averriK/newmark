Dn_BT07 <- function( ky,Sa,Tn,  Mw) {
  # Bray Travasarou 2007 ----
  # Sa in [g]
  # Dn: [cm]
  # Ts <- 1.5 * Tn
  # Ts is not used in this model
  muLnSa <- log(Sa) #Sa(Ts)?  
  
  muLnD <- -1.10 - 2.83 * log(ky) - 0.333 * log(ky)^2 + 0.566 * log(ky) * muLnSa + 3.04 * muLnSa - 0.244 * muLnSa^2 + 0.278 * (Mw - 7)
  sdLnD <- 0.66
  return( data.table(muLnD=muLnD,sdLnD=sdLnD,ID="BT07"))
}
