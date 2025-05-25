Dn_JB07 <- function( PGA, AI=NULL,ky) {
  # Jibson 2007 ----
  # PGA in [g]
  # Dn: [cm]
  # AI in [m/s] >? Confirmar
  # Verificar que AI [m/s] reporta Desplazamientos Dn en [cm] con estos ceficientes.
  # Chequear unidades de AI en paper de Jibson
  if(is.null(AI) & !is.null(PGA) & PGA>0){
    b0 <- 2.6109
    b1 <- 1.9228
    b2 <- 0.0
    AI <- (PGA^b1) * exp(b0 + (b2 / PGA)) # [m/s]
  }
  r <- (ky/PGA)

  muLog10D <- 0.561 * log10(AI) - 3.833 * log10(r) - 1.474
  sdLog10D <- 0.616
  muLnD <- muLog10D*log(10)
  sdLnD <- sdLog10D*log(10)
  return( data.table(muLnD=muLnD,sdLnD=sdLnD,ID="JB07"))
}
