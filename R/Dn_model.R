Dn_model <- function(PGA, ky=NULL, Tn=NULL, Sa=NULL,Mw=NULL, PGV=NULL,AI=NULL, kymin = 0.005, kymax = 0.5, n = 30){
  
  
  if(is.null(ky)){
    ky <- seq(from = log(kymin), to = log(kymax), length.out = n) |>
      exp() |>
      round(digits = 4) |>
      unique()
  }
  
  DT[,ky:=ky]
  DT[,Dm:=exp(muLnD+1/2*((sdLnD)^2))]
  
}
