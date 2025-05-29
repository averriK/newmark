# nolint start
rm(list = ls())
library(data.table)
source("R/Dn_model.R")


Vs30_TARGET <- 300
TR_TARGET <- 2475

UHSTable <- readRDS("data/UHSTable.Rds")
#test
UHS <- UHSTable[TR==TR_TARGET &   p=="mean" & Vref == 760 & Vs30 ==Vs30_TARGET,.(Sa,Tn)]
Ts <- 0.103
ky <- 0.23
Sa <- UHS$Sa
Tn <- UHS$Tn
PGA <- UHS[Tn == 0]$Sa
Dn_model(Sa=Sa,Tn=Tn,PGA=PGA, ky=ky, Ts=Ts, Mw=7.5)


# Robust estimate
# UHS <- UHSTable[TR==TR_TARGET &    Vref == 760 & Vs30 ==Vs30_TARGET,.(Sa,Tn,p)]

