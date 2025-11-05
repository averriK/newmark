# rm(list=ls()) # NO. it deletes params$.
Mw_TARGET <- 6.5
Rrup_gmdp <- 100  # distancia de ruptura efectiva (puede ajustarse)

Vs30_gmdp <- NULL # No additional site classes
#
NS <- 200 # Number of Sampling
# Geometrias. Slopes
Hs <- seq(50, 130, length.out = 5)
s.min <- 2.5
s.max <- 3.5
lo.max <- 0.075 # 0.075
lo.min <- 0.35 # 0.05

# Scenarios
uscs <- list()
uscs[["U1"]] <- c("ML", "MH") # Fine Tailings # Dry Stack
uscs[["U2"]] <- c("SP", "SM") # Coarse Tailings
uscs[["U3"]] <- c("SP", "SM", "ML", "MH") # Coarse Tailings

uscs[["U4"]] <- c("GP", "GW") # Rockfill
uscs[["U5"]] <- c("GP", "GW", "SP", "SM") # Rockfill
uscs[["U6"]] <- c("GP", "GW", "ML", "MH") # Rockfill

uscs[["U7"]] <- c("CL", "CH") # Fine Tailings # Dry Stack
uscs[["U8"]] <- c("CL", "CH", "ML", "MH") # Fine Tailings # Dry Stack
IDm <- paste0("U", seq(1, length(uscs)))

# Newmark Displacement Units
Dn_units <- "mm"
# Target Newmark Displacement (for kmax/kh) in Dn_units
Da_gmdp <- c(10, 100, 1000) # [mm]
# Service Levels
Vref_gmdp <- 760
TR_gmdp <- c(475, 975, 1975, 2475, 4975, 9975)
ID_max <- "max" #
ID_gmdp <- c("gem", "sdLnSaFC1")

# L LOG SCALE
Sa_LOG <- FALSE
Tn_LOG <- FALSE
TR_LOG <- TRUE
