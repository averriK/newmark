rm(list=ls())
library(sdQ)
library(data.table)
library(Hmisc)
source("R/Dn_model.R")
ky <- 0.23
Ts <- 0.103
Mw <- 7.8
Vs30_TARGET <- 300
TR_TARGET <- 2475

UHSTable <- readRDS("data/UHSTable.Rds")
UHS <- UHSTable[TR==TR_TARGET &   Vref == 760 & Vs30 ==Vs30_TARGET,.(Sa,Tn,p)]
PGATable <- UHS[Tn == 0]
UHS[Tn==0, Tn := 0.01,by=.(p)]
SaTable <- UHS[,.(Sa=stats::approx(x = log(Tn), y = log(Sa), xout = log(1.3*Ts))$y |> exp()),by=.(p)]
# 
muLnSa <- log(SaTable[p=="mean"]$Sa)
p <-  as.numeric(SaTable[p!="mean"]$p)
q <-  log(SaTable[p!="mean"]$Sa)
sdLnSa <- sdQ(
  meanValue = muLnSa,
  p = p,
  q = q
)
# Sample Sa

Sa <- rnormQ(
  meanValue = muLnSa,
  p = p,
  q = q,
  n = 1000
) |> exp()

# 
muLnPGA <- log(PGATable[p=="mean"]$Sa)
p <-  as.numeric(PGATable[p!="mean"]$p)
q <-  log(PGATable[p!="mean"]$Sa)
sdLnPGA <- sdQ(
  meanValue = muLnPGA,
  p = p,
  q = q
)
# Sample Sa

PGA <- rnormQ(
  meanValue = muLnPGA,
  p = p,
  q = q,
  n = length(Sa)
) |> exp()

PGV <- (PGA^1.0529) * exp(0.1241) * 100      # PGV back-calculated [cm/s]
AI <- (PGA^1.9228) * exp(2.6109)           # AI back-calculated [m/s]

n <-  length(Sa)

# PGA-based models
# AM88
AUX <- Dn_AM88(PGA = PGA, ky = ky)
AM88Table <- AUX[,.(ID,LnD=rnorm(n=n,mean=muLnD,sd=sdLnD)),by=.(.I)][,.(sample=.I,ID,Dn=exp(LnD))]

#
AUX <- Dn_YG91(PGA = PGA, ky = ky)
YG91Table <- AUX[,.(ID,LnD=rnorm(n=n,mean=muLnD,sd=sdLnD)),by=.(.I)][,.(sample=.I,ID,Dn=exp(LnD))]


# PGA&IA-based models. ky mut be scalar
AUX <- Dn_JB07(PGA = PGA, AI=AI,ky = ky)
JB07Table <- AUX[,.(ID,LnD=rnorm(n=n,mean=muLnD,sd=sdLnD)),by=.(.I)][,.(sample=.I,ID,Dn=exp(LnD))]

#
AUX <- Dn_SR08(PGA = PGA, AI=AI,ky = ky)
SR08Table <- AUX[,.(ID,LnD=rnorm(n=n,mean=muLnD,sd=sdLnD)),by=.(.I)][,.(sample=.I,ID,Dn=exp(LnD))]
# 
AUX <- Dn_BT07(Ts = Ts, Sa = Sa, Mw = Mw, ky = ky)
BT07Table <- AUX[,.(ID,LnD=rnorm(n=n,mean=muLnD,sd=sdLnD)),by=.(.I)][,.(sample=.I,ID,Dn=exp(LnD))]

# BM17
AUX <- Dn_BM17(Ts = Ts, Sa = Sa, Mw = Mw, ky = ky)
BM17Table <- AUX[,.(ID,LnD=rnorm(n=n,mean=muLnD,sd=sdLnD)),by=.(.I)][,.(sample=.I,ID,Dn=exp(LnD))]

#
AUX <- Dn_BM19(Ts = Ts, Sa = Sa, PGA=PGA,PGV=PGV,Mw = Mw, ky = ky)
BM19Table <- AUX[,.(ID,LnD=rnorm(n=n,mean=muLnD,sd=sdLnD)),by=.(.I)][,.(sample=.I,ID,Dn=exp(LnD))]



AUX <- list(AM88Table,YG91Table,JB07Table,SR08Table,BT07Table,BM17Table,BM19Table) |>
  rbindlist(use.names = TRUE, fill = TRUE)
setkey(AUX, ID)                  # fast row access by ID


# LnD <- dcast(AUX, sample ~ ID, value.var = "LnD")

models <- c("YG91","AM88","JB07","BT07","SR08","BM17","BM19")
score <- c(1,3,3,4,4,5,5)
weights <- data.table(
  ID = models,
  weight = score / sum(score)
)
AUX <- weights[AUX, on="ID"]

# Weighted quantiles
# Build Q/uantiles
probs <- c(0.02, 0.05, 0.10, 0.16, 0.50, 0.84, 0.90, 0.95, 0.98)
# DnTable <- AUX[,.(p = as.character(probs),
#                   Dn = quantile(Dn, probs = probs, na.rm = TRUE, type = 7))  # exp after quantile → avoids exponentiating 1 M numbers
# ]
# list(DnTable,data.table(p="mean",Dn=mean(AUX$Dn))) |> rbindlist()

DnTable <- AUX[, .(
  p = probs,
  Dn = wtd.quantile(x= Dn,
    weights = weight,
    probs   = probs,
    type    = 'quantile',
    na.rm   = TRUE
  )
)]
DnTable <- list(DnTable,data.table(p="mean",Dn=mean(AUX$Dn))) |> rbindlist()



# Approach 2: Epistemic and Aleatory uncertainty
