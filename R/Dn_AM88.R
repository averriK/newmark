Dn_AM88 <- function( PGA,  ky) {
  # Ambraseys Menu 1988 ----
  # PGA in [g]
  # Dn: [cm]
  
  r <- (ky/PGA)
  a1 <- 0.90
  a2 <- 2.53
  a3 <- -1.09
  b <- (1 - r)^a2
  c <- r^a3
  
  # MODEL ----
  muLog10D <- a1 + log10(b * c)
  sdLog10D <- 0.30
  muLnD <- muLog10D*log(10)
  sdLnD <- sdLog10D*log(10)
  return( data.table(muLnD=muLnD,sdLnD=sdLnD,ID="AM88"))
  
}
