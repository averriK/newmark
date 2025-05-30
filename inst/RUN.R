
rm(list=ls())
library(data.table)
# library(newmark)     # ← your package
devtools::load_all()  # load package in development mode
## ---- LOAD & FILTER UHS ----------------------------------------------------
Vs30_TARGET <- c(300, 475)    # ← your target array
TR_TARGET <- c(2475,475)
## ---- 1  Load UHSTable (no filtering on TR now)  -------------
UHSTable <- readRDS("inst/UHSTable.Rds")

UHS <- UHSTable[TR %in%  TR_TARGET&  Vref == 760 &    Vs30 == 760,              # reference-rock rows only
  .(TR, Sa, Tn, p,Vref)
]

## ---- 2  SaF for every (TR, Vs30) ----------------------------
SaF_all <- fitSaF(
  uhs  = UHS,
  vs30 = Vs30_TARGET,
  NS   = 100,
  vref = 760
)

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



