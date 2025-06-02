
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
  ][, fitDn(.SD, ky = kyi, Ts = Ts, Mw = 6.5, NS = 1000,uncertainty="none"), by = .(TR, Vs30)][, ky := kyi][Dn < 1e-16, Dn := 0]
  
}) |> rbindlist()
DnTable[TR==475 & Vs30 ==400 & p %in% c(0.16,0.50,"mean",0.84) & ky==0.2]
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
  fitDn(.SD, ky = ky, Ts = Ts, Mw = 6.5, NS = 1000,,uncertainty="Sa"),
  by = .(TR, Vs30)
][abs(Dn) < 1e-16, Dn := 0][]
DnTable[TR==475 & Vs30 ==400 & p %in% c(0.16,0.50,"mean",0.84)]


##

UHS <- UHSTable[ Vref==760 &  Vs30 ==300 & TR==2475 ]
outSa <- sampleSa(UHS, Td=0.10, n=1000)
summary(outSa$Sa)
quantile(outSa$Sa, probs=c(0.01,0.05,0.5,0.95,0.99))
