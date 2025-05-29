
library(sdQ)
library(data.table)
source("R/Dn_model.R")
ky <- 0.23
Ts <- 0.103
Td <- 1.5*Ts
Mw <- 7.8
Vs30_TARGET <- 300
TR_TARGET <- 2475

UHSTable <- readRDS("data/UHSTable.Rds")
UHS <- UHSTable[TR==TR_TARGET &   Vref == 760 & Vs30 ==Vs30_TARGET,.(Sa,Tn,p)]
PGATable <- UHS[Tn == 0]
UHS[Tn==0, Tn := 0.01,by=.(p)]
SaTable <- UHS[,.(Sa=stats::approx(x = log(Tn), y = log(Sa), xout = log(Td))$y |> exp()),by=.(p)]
# 
muLnSa <- log(SaTable[p=="mean"]$Sa)
p <-  as.numeric(SaTable[p!="mean"]$p)
q <-  log(SaTable[p!="mean"]$Sa)
sdLnSa <- sdQ(
  meanValue = muLnSa,
  p = p,
  q = q
)
sdLnSa
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

PGV <- (PGA^1.0529) * exp(0.1241) * 100      # convert to cm / s
AI <- (PGA^1.9228) * exp(2.6109)           # AI back-calculated [m/s]

n <-  length(Sa)

# PGA-based models
# AM88
AUX <- Dn_AM88(PGA = PGA, ky = ky)
AM88Table <- AUX[,.(ID,LnD=rnorm(n=n,mean=muLnD,sd=sdLnD)),by=.(.I)][,.(ID,LnD)]

#
AUX <- Dn_YG91(PGA = PGA, ky = ky)
YG91Table <- AUX[,.(ID,LnD=rnorm(n=n,mean=muLnD,sd=sdLnD)),by=.(.I)][,.(ID,LnD)]


# PGA&IA-based models. ky mut be scalar
AUX <- Dn_JB07(PGA = PGA, AI=AI,ky = ky)
JB07Table <- AUX[,.(ID,LnD=rnorm(n=n,mean=muLnD,sd=sdLnD)),by=.(.I)][,.(ID,LnD)]

#
AUX <- Dn_SR08(PGA = PGA, AI=AI,ky = ky)
SR08Table <- AUX[,.(ID,LnD=rnorm(n=n,mean=muLnD,sd=sdLnD)),by=.(.I)][,.(ID,LnD)]
# 
AUX <- Dn_BT07(Ts = Ts, Sa = Sa, Mw = Mw, ky = ky)
BT07Table <- AUX[,.(ID,LnD=rnorm(n=n,mean=muLnD,sd=sdLnD)),by=.(.I)][,.(ID,LnD)]

# BM17
AUX <- Dn_BM17(Ts = Ts, Sa = Sa, Mw = Mw, ky = ky)
BM17Table <- AUX[,.(ID,LnD=rnorm(n=n,mean=muLnD,sd=sdLnD)),by=.(.I)][,.(ID,LnD)]

#
AUX <- Dn_BM19(Ts = Ts, Sa = Sa, PGA=PGA,Mw = Mw, ky = ky)
BM17Table <- AUX[,.(ID,LnD=rnorm(n=n,mean=muLnD,sd=sdLnD)),by=.(.I)][,.(ID,LnD)]







