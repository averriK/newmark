###############################################################################
# runDn.R ― Secuencial con trazas y saneo estricto de columnas
###############################################################################

plan(sequential)   # future.apply ― ejecución secuencial para depuración

SHT <- ShearTable[level == "mean", .(IDg, IDm, VS30, mo, Ts)]

# ---------------------------------------------------------------------------
# Configuración base‑rock para la amplificación
# ---------------------------------------------------------------------------
Vref_gmdp <- 760L
UHSRock <- UHSTable[
  ID == "gem" & Vs30 == Vref_gmdp & Vref == Vref_gmdp,
  .(ID, TR, Tn, Sa, p, Vs30, Vref)
]

# ---------------------------------------------------------------------------
# Tablas de salida
# ---------------------------------------------------------------------------
DnTable   <- data.table()
TR_TARGET <- TR_gmdp          # vector de TR definidos por el usuario

# ---------------------------------------------------------------------------
# Bucle principal: geometría × sección
# ---------------------------------------------------------------------------
for (gi in unique(SHT$IDg)) {
  for (mi in unique(SHT$IDm)) {
    
    VSi <- round(SHT[IDg == gi & IDm == mi, VS30]); stopifnot(length(VSi) == 1L)
    TSi <- SHT[IDg == gi & IDm == mi, Ts];          stopifnot(length(TSi) == 1L)
    
    cat(sprintf("> Site response spectra (UHS) for Vs30=%g (section: %s geometry:%s)...\n",
                VSi, mi, gi))
    
    # ----------------- Amplificación de sitio (fitSaF) ------------------
    UHS <- UHSRock[
      ,
      newmark::fitSaF(
        .SD[, .(Tn, p, Sa)],   # ⬅︎ SOLO columnas canónicas
        vs30 = VSi,
        vref = Vref_gmdp,
        ns   = NS,
        Rrup = Rrup_gmdp
      )[, .(
        TR   = .BY$TR,
        Vs30 = VSi,
        Vref = Vref_gmdp,
        ID   = "sdLnSaFC1",
        Tn, p, Sa = SaF
      )],
      by = .(TR)
    ]
    
    # --------------------- Bucles TR → ky → fitDn ----------------------
    for (TRi in TR_TARGET) {
      cat(sprintf("> Newmark displacements (Dn) for Ts = %g s  and TR = %d yr...\n",
                  TSi, TRi))
      
      Sa <- UHS[TR == TRi & Tn == min(Tn) & p == "mean"]$Sa
      if (length(Sa) == 0L ) next

      ky_min <- Sa * 0.25
      ky_max <- Sa * 1.25
      if (ky_min >= ky_max) next
      
      ky_TARGET <- seq(ky_min, ky_max, length.out = 20L)
      
      for (k in ky_TARGET) {
        cat(sprintf("  >> Evaluating ky = %.4f\n", k))
        
        out <- tryCatch({
          tmp <- UHS[TR == TRi,
                     newmark::fitDn(
                       .SD[, .(Tn, p, Sa)],   # ⬅︎ columnas limpias
                       ky   = k,
                       Ts   = TSi,
                       Mw   = Mw_TARGET,
                       NS   = NS,
                       Rrup = Rrup_gmdp
                     )]
          tmp[, `:=`(
            TR   = TRi,
            ky   = k,
            Ts   = TSi,
            Vs30 = VSi,
            IDg  = gi,
            IDm  = mi
          )]
        }, error = function(e) {
          cat(sprintf("  !! Error for ky = %.4f: %s\n", k, e$message))
          NULL
        })
        
        if (!is.null(out)) {
          DnTable <- rbindlist(list(DnTable, out), use.names = TRUE)
        }
      }
    }
  }
}

# ---------------------------------------------------------------------------
# Conversión de unidades
# ---------------------------------------------------------------------------
if (Dn_units == "cm") DnTable[, units := "cm"]
if (Dn_units == "mm") { DnTable[, units := "mm"]; DnTable[, Dn := Dn * 10] }
if (Dn_units == "m")  { DnTable[, units := "m"];  DnTable[, Dn := Dn / 100] }
if (Dn_units == "in") { DnTable[, units := "in"]; DnTable[, Dn := Dn / 2.54] }

plan(sequential)   # restaurar plan por si se usa fuera de este script
