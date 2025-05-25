Dn_SR08 <- function(PGA, ky,AI=NULL) {
  # Saygili Rathje 2008
  # PGA in [g]
  # Dn: [cm]
  # AI in [m/s] >? Confirmar
  if(is.null(AI)){
    b0 <- 2.6109
    b1 <- 1.9228
    b2 <- 0.0
    
    AI <- (PGA^b1) * exp(b0 + (b2 / PGA)) # [m/s]
  }
  r <- (ky/PGA)
  
  muLnD <- 2.39 - 5.24 * r - 18.78 * r^2 + 42.01 * r^3 - 29.15 * r^4 - 1.56 * log(PGA) + 1.38 * log(AI)
  sdLnD <- 0.46 + 0.56 * r
  return( data.table(muLnD=muLnD,sdLnD=sdLnD,ID="SR08"))
}
