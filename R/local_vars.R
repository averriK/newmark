# Work around R CMD check “no visible binding” notes
utils::globalVariables(
    c("..ID","AF", "w","pga", "PGA", "N", "depth", "lat", "lon", "SaRock", "TnVal", "FValue", "SaF", "LnF", "muLnF", "sdLnF", "pnum", "lnSa", "Tn", "p", "Dn", "weight", "muLnD", "sdLnD", "LnD", "LnSaF", "muLnSaF", "sdLnSaF", "ID", "Sa", "SaF_SS17", "TR", "Vs30", "Vref", ".", "AEP", "IT", "POE", "iml", "imt")
)
