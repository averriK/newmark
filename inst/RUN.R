
rm(list=ls())
library(data.table)
# library(newmark)     # ← your package
devtools::load_all()  # load package in development mode
## ---- LOAD & FILTER UHS ----------------------------------------------------
Vs30_TARGET <- 300    # ← your target array
TR_TARGET <- 2475
## ---- 1  Load UHSTable (no filtering on TR now)  -------------
UHSTable <- readRDS("inst/UHSTable.Rds")
# 
# UHS <- UHSTable[TR %in%  TR_TARGET&  Vref == 760 & Vs30 %in% Vs30_TARGET& !(p %in% c("0.95","0.98","0.05","0.02")),              # reference-rock rows only
#   .(TR, Sa, Tn, p,Vref,Vs30)
# ]

# & !(p %in% c("0.95","0.98","0.05","0.02"))


Ts <- 0.13   # fundamental period of sliding mass (s)
ky_TARGET <- seq(0.05, 0.5, by = 0.025)

DnTable <- lapply(ky_TARGET, function(kyi) {
  # For each kyi, call the same logic as before:
 UHSTable[
    TR %in% c(475, 2475) &
      Vs30 %in% c(300, 400) &
      Vref == 760
  ][, fitDn(.SD, ky = kyi, Ts = Ts, Mw = 6.5, NS = 100,uncertainty="none"), by = .(TR, Vs30)][, ky := kyi][Dn < 1e-16, Dn := 0]
  
}) |> rbindlist()

## ------------------------------------------------------------
## 2)  Build Newmark-displacement quantiles
##     (fitDn uses the same UHS; no extra pre-processing needed)
## ------------------------------------------------------------
Ts <- 0.13   # fundamental period of sliding mass (s)
ky<- 0.15  # yield acceleration (g)

DnTable <- UHSTable[
  TR %in% c(475, 2475) &
    Vs30 %in% c(300, 400) &
    Vref == 760
][,
  fitDn(.SD, ky = ky, Ts = Ts, Mw = 6.5, NS = 1000),
  by = .(TR, Vs30)
][abs(Dn) < 1e-16, Dn := 0][]

