
rm(list=ls())
library(data.table)
# library(newmark)     # ← your package
devtools::load_all()  # load package in development mode
## ---- LOAD & FILTER UHS ----------------------------------------------------
UHSTable <- readRDS("inst/UHSTable.Rds")

UHS <- UHSTable[
    Vref == 760 &
    Vs30 == 760,
  .(TR,Sa, Tn, p)   # keep only what fitDn() needs
]
## ---- 2  SaF quantiles for every TR  ------------------------------------

SaF_all <- rbindlist(
  lapply(unique(UHS$TR), function(tr) {
    
    df <- fitSaF(
      uhs  = UHS[TR == tr, .(Sa, Tn, p)],   # pass only the columns fitSaF needs
      vs30 = 800,
      NS   = 1000,
      vref = 760
    )
    
    df[, TR := tr]                     # tag result with current TR
    df
  }),
  use.names = TRUE
)

setcolorder(SaF_all, c("TR", "Tn", "p", "SaF"))

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



