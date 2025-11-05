# runKmax.R
p_TARGET <- "mean" # confidence level ("mean" or numeric)
Tn_PGA <- UHSTable$Tn |> min()
# PGA table for rock reference --------------------------------------------
Vref_gmdp <- 760
PGATable <- UHSTable[
    ID == "gem" & Tn == Tn_PGA & Vs30 == Vref_gmdp & Vref == Vref_gmdp,
    .(PGA = SaF),
    by = .(TR, p)
]

# merge PGA into DnTable ---------------------------------------------------
setkey(DnTable, TR, p)
AUX <- DnTable[PGATable, on = .(TR, p), allow.cartesian = TRUE]

# kh coefficients ----------------------------------------------------------

# kmax coefficients --------------------------------------------------------


kmaxTable <- DnTable[
  , kmaxMC(Dn, ky, p, Da = Da_gmdp, p_TARGET = p_TARGET, ns = 1000),
  by = .(TR, IDg, IDm)
][, p := p_TARGET]
