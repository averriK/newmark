# ============================================================
#  Public wrapper: fitSaF()
#  (loops over TR, vs30, and all periods in each UHS slice)
# ============================================================
# nolint start
#' Quantiles of *SaF* = Sa × F(ST17) for:
#'   • every return period in the input UHSTable
#'   • every Vs30 value supplied (scalar or vector)
#'   • every oscillator period stored in the UHS
#'
#' @param uhs   data.table – FULL UHSTable (must contain **TR**, **Sa**, **Tn**, **p**).
#' @param vs30  numeric scalar **or vector** – target Vs30 values (m/s).
#' @param NS    integer ≥ 1 – Monte-Carlo samples per (TR, Vs30, Tn) combo (default 30).
#' @param vref  numeric scalar – reference Vs (default 760 m/s).
#'
#' @return data.table with columns **TR**, **Vs30**, **Tn**, **p**, **SaF**.
#' @export

fitSaF <- function(uhs,
                   vs30,
                   NS   = 30,
                   vref = 760) {
  
  stopifnot(all(c("TR", "Sa", "Tn", "p") %in% names(uhs)))
  
  TR_vals <- sort(unique(uhs$TR))
  VS_vals <- sort(unique(vs30))
  
  OUT <- rbindlist(
    lapply(VS_vals, function(VS) {
      
      rbindlist(
        lapply(TR_vals, function(TR_i) {
          
          UHS_TR  <- uhs[TR == TR_i, .(Sa, Tn, p)]
          UHS_mod <- copy(UHS_TR)
          UHS_mod[Tn == 0, Tn := 0.01, by = .(p)]
          
          RES_TR <- rbindlist(
            lapply(sort(unique(UHS_TR$Tn)), function(Tn_orig) {
              
              Tn_eval <- if (Tn_orig == 0) 0.01 else Tn_orig
              Sa_Tn   <- sampleSa(UHS_mod, Td = Tn_eval, n = NS)$Sa
              pga     <- sampleSa(UHS_mod, Td = 0.01,   n = NS)$Sa
              
              SaF_tbl <- getSaF(
                SaF_ST17,
                Sa   = Sa_Tn,
                pga  = pga,
                Tn   = Tn_orig,
                vs30 = VS,
                vref = vref,
                n    = NS)
              
              probs  <- UHS_TR[p != "mean", unique(as.numeric(p))]
              SaF_q  <- data.table(
                p   = probs,
                SaF = stats::quantile(SaF_tbl$SaF,
                                      probs = probs,
                                      names = FALSE,
                                      type  = 7))
              SaF_mu <- data.table(p = "mean", SaF = mean(SaF_tbl$SaF))
              
              rbindlist(list(SaF_q, SaF_mu))[, Tn := Tn_orig]
            }),
            use.names = TRUE)
          
          RES_TR[, `:=`(TR = TR_i, Vs30 = VS)]
          RES_TR
        }),
        use.names = TRUE)
    }),
    use.names = TRUE)
  
  OUT[, Vref := vref]
  setcolorder(OUT, c("TR", "Vs30", "Vref", "Tn", "p", "SaF"))
  OUT[]
}



# nolint end
