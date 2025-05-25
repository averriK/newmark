rm(list=ls())
library(data.table)
source("R/Dn_AM88.R")
source("R/Dn_BM17.R")
source("R/Dn_BM19.R")
source("R/Dn_BT07.R")
source("R/Dn_JB07.R")
source("R/Dn_SR08.R")
source("R/Dn_YG91.R")




Dn_AM88(PGA=0.7,ky=0.23)
Dn_BT07(Sa=0.931, Mw=8.5,ky=0.23)
Dn_BM17(Tn=0.06, Sa=0.931, Mw=8.5,ky=0.23) # Subduction
Dn_BM19(PGA=0.7, Tn=0.06, Sa=0.931, Mw=8.5,ky=0.23) # Shallow Crustal
Dn_JB07(PGA=0.7,ky=0.23)
Dn_SR08(PGA=0.7,ky=0.23)
Dn_YG91(PGA=0.7,ky=0.23)
