
rm(list=ls())
library(data.table)
# library(newmark)     # ← your package
devtools::load_all()  # load package in development mode
## ---- LOAD & FILTER UHS ----------------------------------------------------
UHSTable <- readRDS("inst/UHSTable.Rds")

UHS <- UHSTable[
  TR   == 2475 &
    Vref == 760 &
    Vs30 == 760,
  .(Sa, Tn, p)   # keep only what fitDn() needs
]

## ------------------------------------------------------------
## 1)  Build SaF(Tn) quantiles for a single Vs30
## ------------------------------------------------------------

SaF <- fitSaF(
  uhs  = UHS,
  vs30 = 800,
  Tn   = 0.50,
  NS   = 1000,
  vref = 760
)

lapply(UHS$Tn |> unique(), function(Tn) {
  fitSaF(
    uhs  = UHS,
    vs30 = 800,
    Tn   = Tn,
    NS   = 1000,
    vref = 760
  )
}) |> rbindlist()
## ------------------------------------------------------------
## 2)  Build Newmark-displacement quantiles
##     (fitDn uses the same UHS; no extra pre-processing needed)
## ------------------------------------------------------------
Ts <- 1.0   # fundamental period of sliding mass (s)
ky <- 0.15  # yield acceleration (g)

Dn <- fitDn(
  uhs = UHS,
  ky  = 0.15,
  Ts  = 1.0,
  Mw  = 6.5,
  NS  = 1000
)



