# runASCE.R
## -----------------------------------------------------------------
##  build MCER spectra (TR = 2475) for both “gem” and “gmdp”
## -----------------------------------------------------------------

DT1 <- UHSTable[ID == ID_max & p == "mean",
                {
                  newmark::designUHS(UHS = .SD, TL = 8,spectrum = "MCER")
                },
                by = .(Vs30, ID, p, TR) 
][, ID := "mcer"]

DT2 <- UHSTable[ID == ID_max & p == "mean",
         {
           newmark::designUHS(UHS = .SD, TL = 8,spectrum = "DESIGN")
         },
         by = .(Vs30, ID, p, TR) 
][, ID := "design"]
ASCETable <- list(DT1,DT2) |> rbindlist(use.names = TRUE)

setorder(ASCETable, TR, Vs30, p, ID, Tn)
 C
