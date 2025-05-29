# nolint start
rm(list = ls())

library(sdQ)
library(data.table)
library(Hmisc)

source("R/Dn_model.R")
source("R/helpers.R")

# Data read from file. Should be an ARGUMENT from my funciton (UHSTable)
UHSTable <- readRDS("data/UHSTable.Rds")
UHS <- UHSTable[TR == TR_TARGET & Vref == 760 & Vs30 == Vs30_TARGET, .(Sa, Tn, p)]
ky <- 0.16
Ts <- 0.103
Mw <- 7.8
NS <- 2000
Vs30_TARGET <- 300
TR_TARGET <- 2475
BM_model <- "crustal" # c("shallow,"crustal","interface","subduction")
