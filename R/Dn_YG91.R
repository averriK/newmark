Dn_YG91 <- function( PGA, ky) {
  # Yegian 1991 ----
  # PGA in [g]
  # Dn: [cm]
  r <- (ky/PGA)
  muLog10D <- 0.22 - 10.12 * r + 16.38 * r^2 - 11.48 * r^3
  sdLog10D <- 0.45
  muLnD <- muLog10D*log(10)
  sdLnD <- sdLog10D*log(10)
  
  return( data.table(muLnD=muLnD,sdLnD=sdLnD,ID="YG91"))
}
