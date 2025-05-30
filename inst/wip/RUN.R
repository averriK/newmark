#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## RUN.R – quick functional test for fitDn()
## ---------------------------------------------------------------------------
# devtools::load_all() # load package in development mode  
suppressPackageStartupMessages({
  library(data.table)          # data wrangling
library(sdQ)
})
source("R/api.R")
source("R/helpers.R")

## ---- INPUT SETTINGS -------------------------------------------------------
Vs30_TARGET <- 300     # m/s
TR_TARGET   <- 2475    # return period (years)
BM_model    <- "crustal"   # "crustal" / "interface" branch for Bray & Macedo
ky          <- 0.16        # yield acceleration (g)
Ts          <- 0.103       # fundamental period of sliding mass (s)
Mw          <- 7.8         # scenario magnitude
NS          <- 1000        # Monte-Carlo samples per model
set.seed(1234)             # reproducibility

## ---- LOAD & FILTER UHS ----------------------------------------------------
UHSTable <- readRDS("data/UHSTable.Rds")

UHS <- UHSTable[
  TR   == TR_TARGET &
    Vref == 760 &
    Vs30 == Vs30_TARGET,
  .(Sa, Tn, p)   # keep only what fitDn() needs
]

if (nrow(UHS) == 0L)
  stop("No rows matched TR_TARGET / Vs30_TARGET – check inputs.")

## ---- RUN DISPLACEMENT ENSEMBLE -------------------------------------------
fitDn(
  UHS      = copy(UHS),     # copy() so original isn’t modified in place
  ky       = ky,
  Ts       = Ts,
  Mw       = Mw,
  NS       = NS,
  BM_model = BM_model
)
