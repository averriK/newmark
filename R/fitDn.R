
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

n <-  length(Sa)
muY <- double(n)
sdY <- double(n)
i <- 1
for (i in seq_len(n)) {
  OUT <- Dn_model(Ts = Ts, Sa = Sa[i], PGA=PGA[i],Mw = Mw, ky = ky)
  
  muY[i] <- OUT$muLnD
  sdY[i] <- OUT$sdLnD
  # ID[i] <- OUT$ID
}











