# nolint start
rm(list = ls())

library(sdQ)
library(data.table)
library(Hmisc)

source("R/Dn_models.R")
source("R/helpers.R")
ky <- 0.16
Ts <- 0.103
Mw <- 7.8
NS <- 2000
Vs30_TARGET <- 300
TR_TARGET <- 2475
models <- c("YG91", "AM88", "JB07", "BT07", "SR08", "BM17", "BM19")
score <- c(1, 2, 2, 3, 3, 4, 4)
weights <- data.table(
    ID = models,
    weight = score / sum(score)
)

BM_model <- "crustal" # c("shallow,"crustal","interface","subduction")


# Data read from file. Should be an ARGUMENT from my funciton (UHSTable)
UHSTable <- readRDS("data/UHSTable.Rds")


# Step 0. UHS Table ----
UHS <- UHSTable[TR == TR_TARGET & Vref == 760 & Vs30 == Vs30_TARGET, .(Sa, Tn, p)]
UHS[Tn == 0, Tn := 0.01, by = .(p)]





# Step 1. PGA Uncertainty -----

# PGATable <- UHS[Tn == 0.01]
# muLnPGA <- log(PGATable[p == "mean"]$Sa)
# p <- as.numeric(PGATable[p != "mean"]$p)
# q <- log(PGATable[p != "mean"]$Sa)
# sdLnPGA <- sdQ(
#     meanValue = muLnPGA,
#     p = p,
#     q = q
# )
# PGA <- rnormQ(
#     meanValue = muLnPGA,
#     p = p,
#     q = q,
#     n = NS
# ) |> exp()

OUT <- sample_Sa(UHS, Td = 0.01, n = NS)
PGA <- OUT$Sa

PGV <- (PGA^1.0529) * exp(0.1241) * 100 # PGV back-calculated [cm/s]
AI <- (PGA^1.9228) * exp(2.6109) # AI back-calculated [m/s]


# Step 2: Sa Uncertainty for Td=Ts ----

# Td <- 1.0 * Ts
# SaTable <- UHS[, .(Sa = stats::approx(x = log(Tn), y = log(Sa), xout = log(Td))$y |> exp()), by = .(p)]
# 
# muLnSa <- log(SaTable[p == "mean"]$Sa)
# p <- as.numeric(SaTable[p != "mean"]$p)
# q <- log(SaTable[p != "mean"]$Sa)
# sdLnSa <- sdQ(
#     meanValue = muLnSa,
#     p = p,
#     q = q
# )
# Sa_10_Ts <- rnormQ(
#     meanValue = muLnSa,
#     p = p,
#     q = q,
#     n = NS
# ) |> exp()
OUT <- sample_Sa(UHS, Td = 1.0 * Ts, n = NS)
Sa_10_Ts <- OUT$Sa


# Step 3: Sa Uncertainty for Td=1.3*Ts ----

# Td <- 1.3 * Ts
# SaTable <- UHS[, .(Sa = stats::approx(x = log(Tn), y = log(Sa), xout = log(Td))$y |> exp()), by = .(p)]
# 
# muLnSa <- log(SaTable[p == "mean"]$Sa)
# p <- as.numeric(SaTable[p != "mean"]$p)
# q <- log(SaTable[p != "mean"]$Sa)
# sdLnSa <- sdQ(
#     meanValue = muLnSa,
#     p = p,
#     q = q
# )
# Sa_13_Ts <- rnormQ(
#     meanValue = muLnSa,
#     p = p,
#     q = q,
#     n = NS
# ) |> exp()
OUT <- sample_Sa(UHS, Td = 1.3*Ts, n = NS)
Sa_13_Ts <- OUT$Sa


# Step 4: Sa Uncertainty for Td=1.5*Ts ----

Td <- 1.5 * Ts
SaTable <- UHS[, .(Sa = stats::approx(x = log(Tn), y = log(Sa), xout = log(Td))$y |> exp()), by = .(p)]

muLnSa <- log(SaTable[p == "mean"]$Sa)
p <- as.numeric(SaTable[p != "mean"]$p)
q <- log(SaTable[p != "mean"]$Sa)
sdLnSa <- sdQ(
    meanValue = muLnSa,
    p = p,
    q = q
)
Sa_15_Ts <- rnormQ(
    meanValue = muLnSa,
    p = p,
    q = q,
    n = NS
) |> exp()
# Sa_15_Ts <- OUT$Sa



# Step 5:  Newmark Displacements for PGA-Based models   ----
DnTable <- data.table()
# AM88
# AUX <- Dn_AM88(PGA = PGA, ky = ky)
DnSample <- Dn_model(.fun = Dn_AM88, PGA = PGA, ky = ky, n = NS)
# DnSample <- AUX[,.(ID,LnD=rnorm(n=NS,mean=muLnD,sd=sdLnD)),by=.(.I)][,.(sample=.I,ID,Dn=exp(LnD))]
DnTable <- list(DnTable, DnSample) |> rbindlist(use.names = TRUE, fill = TRUE)
#
# AUX <- Dn_YG91(PGA = PGA, ky = ky)
DnSample <- Dn_model(.fun = Dn_YG91, PGA = PGA, ky = ky, n = NS)
# DnSample <- AUX[,.(ID,LnD=rnorm(n=NS,mean=muLnD,sd=sdLnD)),by=.(.I)][,.(sample=.I,ID,Dn=exp(LnD))]
DnTable <- list(DnTable, DnSample) |> rbindlist(use.names = TRUE, fill = TRUE)

# PGA&IA-based models. ky mut be scalar
# AUX <- Dn_JB07(PGA = PGA, AI=AI,ky = ky)
DnSample <- Dn_model(.fun = Dn_JB07, PGA = PGA, AI = AI, ky = ky, n = NS)
# DnSample <- AUX[,.(ID,LnD=rnorm(n=NS,mean=muLnD,sd=sdLnD)),by=.(.I)][,.(sample=.I,ID,Dn=exp(LnD))]
DnTable <- list(DnTable, DnSample) |> rbindlist(use.names = TRUE, fill = TRUE)

#
# AUX <- Dn_SR08(PGA = PGA, AI=AI,ky = ky)
DnSample <- Dn_model(.fun = Dn_SR08, PGA = PGA, AI = AI, ky = ky, n = NS)
# DnSample <- AUX[,.(ID,LnD=rnorm(n=NS,mean=muLnD,sd=sdLnD)),by=.(.I)][,.(sample=.I,ID,Dn=exp(LnD))]
DnTable <- list(DnTable, DnSample) |> rbindlist(use.names = TRUE, fill = TRUE)

# Step 6: Newmark Displacements for Sa-Based models   ----
# AUX <- Dn_BT07(Ts = Ts, Sa = Sa_10_Ts, Mw = Mw, ky = ky)
DnSample <- Dn_model(.fun = Dn_BT07, Ts = Ts, Sa = Sa_10_Ts, Mw = Mw, ky = ky, n = NS)
# DnSample <- AUX[,.(ID,LnD=rnorm(n=NS,mean=muLnD,sd=sdLnD)),by=.(.I)][,.(sample=.I,ID,Dn=exp(LnD))]
DnTable <- list(DnTable, DnSample) |> rbindlist(use.names = TRUE, fill = TRUE)

# BM17

# BM_model <- "crustal" # c("shallow,"crustal","interface","subduction")

if (tolower(BM_model) %in% c("subduction", "interface")) {
    # Subduction
    # AUX <- Dn_BM17(Ts = Ts, Sa = Sa_15_Ts, Mw = Mw, ky = ky)
    DnSample <- Dn_model(.fun = Dn_BM17, Ts = Ts, Sa = Sa_15_Ts, Mw = Mw, ky = ky, n = NS)
    # DnSample <- AUX[,.(ID,LnD=rnorm(n=NS,mean
    DnTable <- list(DnTable, DnSample) |> rbindlist(use.names = TRUE, fill = TRUE)
}

if (tolower(BM_model) %in% c("shallow", "crustal")) {
    # Crustal
    # AUX <- Dn_BM19(Ts = Ts, Sa = Sa_13_Ts, PGA=PGA,PGV=PGV,Mw = Mw, ky = ky)
    DnSample <- Dn_model(.fun = Dn_BM19, Ts = Ts, Sa = Sa_13_Ts, PGA = PGA, PGV = PGV, Mw = Mw, ky = ky, n = NS)
    # DnSample <- AUX[,.(ID,LnD=rnorm(n=NS,mean=muLnD,sd=sdLnD)),by=.(.I)][,.(sample=.I,ID,Dn=exp(LnD))]
    DnTable <- list(DnTable, DnSample) |> rbindlist(use.names = TRUE, fill = TRUE)
}
# Step 7: Weighted quantiles ----


AUX <- weights[DnTable, on = "ID"]
AUX[, units := "cm"]

probs <- UHS[p != "mean"]$p |> unique() |> as.numeric()

DnTable <- AUX[, .(
    p = probs,
    Dn = Hmisc::wtd.quantile(
        x = Dn,
        weights = weight,
        probs = probs,
        type = "quantile",
        na.rm = TRUE
    )
)]
DnTable <- list(DnTable, data.table(p = "mean", Dn = mean(AUX$Dn))) |> rbindlist()
