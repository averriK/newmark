rm(list=ls())
library(data.table)
# library(newmark)     # ← your package
devtools::load_all()  # load package in development mode

## ---- LOAD & FILTER UHS ----------------------------------------------------
Vs30_TARGET <- 300    # ← your target array
TR_TARGET <- 2475

## ---- 1  Load UHSTable (no filtering on TR now) -------------
UHSTable <- readRDS("inst/UHSTable.Rds")
Ts <- 0.13   # fundamental period of sliding mass (s)
# We'll show two tests of fitDn() and then one test of fitSaF().

# ---------------------------------------------------------------------------
# 1) EXAMPLE 1: Loop over multiple ky values, "uncertainty = 'none'"
# ---------------------------------------------------------------------------


DnTable_loop <- lapply(seq(0.05, 0.5, by = 0.025), function(kyi) {
  # For each kyi, call fitDn()
  UHSTable[
    TR %in% c(475, 2475) &
      Vs30 %in% c(300, 400) &
      Vref == 760
  ][
    , fitDn(
      .SD, 
      ky   = kyi, 
      Ts   = Ts, 
      Mw   = 6.5, 
      NS   = 1000,
      uncertainty = "none"   # hazard = mean, Dn = mean
    ),
    by = .(TR, Vs30)
  ][
    # clamp near-zero displacements
    Dn < 1e-16, Dn := 0
  ][
    # add a column for the current ky
    , ky := kyi
  ]
}) |> rbindlist()

# Quick peek at a subset:
DnTable_loop[ R==475 & Vs30 ==400 & p %in% c(0.16,0.50,"mean",0.84) & ky==0.2]

# ---------------------------------------------------------------------------
# 2) EXAMPLE 2: Single call with "uncertainty = c('Sa','Dn')"
# ---------------------------------------------------------------------------
ky <- 0.15
DnTable_both <- UHSTable[
  TR %in% c(475, 2475) &
    Vs30 %in% c(300, 400) &
    Vref == 760
][
  , fitDn(
    .SD, 
    ky = ky, 
    Ts = Ts, 
    Mw = 6.5, 
    NS = 1000,
    uncertainty = c("Sa","Dn")  # hazard + displacement randomness
  ),
  by = .(TR, Vs30)
][
  abs(Dn) < 1e-16, 
  Dn := 0
][]

DnTable_both[TR==475 & Vs30 ==400 & p %in% c(0.16,0.50,"mean",0.84)]

# ---------------------------------------------------------------------------
# 3) EXAMPLE 3: Test fitSaF() for site factor
# ---------------------------------------------------------------------------
# We'll pick a single TR, vref=760 => "rock" subset, 
# and call fitSaF with vs30=300, "uncertainty='both'"

uhs_sub <- UHSTable[ TR==475 & Vref==760 & Vs30==760 ]  # single TR, single vref=760
# must contain Tn=0.01 for each p => used as "PGA"
# must contain Tn=0.1,0.2,... etc. for each p => used as "Sa"
# plus p="mean" row.

res <- fitSaF(
  uhs_sub,
  vs30         = 300,
  NS           = 1000,
  models       = "ST17",
  score        = 1,
  uncertainty  = c("Sa","F")
)

res[Tn==0.2 & p %in% c("0.16","0.5","mean","0.84")]
