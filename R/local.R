# Work around R CMD check “no visible binding” notes
utils::globalVariables(
  c("pnum","lnSa","Tn", "p", "Dn", "weight", "muLnD", "sdLnD", "LnD", "LnSaF", "muLnSaF", "sdLnSaF","ID","Sa" ,"SaF_SS17" ,"TR", "Vs30", "Vref","."))
