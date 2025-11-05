# runUHS.R
# =====================================================================
#  Builds UHS & site‑factor tables for every combination of:
#    • TR       : 475–9975
#    • Vs30     : as defined in gem folders and user-defined list
#    • Vref     : 760 or 3000, depending on gem_path content
#    • ID       : model name for rock, "model.tag" for site-amplified
# =====================================================================

plan(multisession, workers = max(1L, parallel::detectCores() - 1L))

UHSTable <- data.table::data.table()
AEPTable <- data.table::data.table()

# ---------------------------------------------------------------------
# 0. Check parameters
# ---------------------------------------------------------------------
if(!exists("Rrup_gmdp")) Rrup_gmdp <- 100
if(!exists("ID_max")) ID_max <- "max"
if(!exists("Vs30_gmdp")) Vs30_gmdp <- c(180,270,360,560,760,1250)  # fixed typo
if(!exists("gem_path")) gem_path <- file.path("examples","oq","gem")
if(!exists("TR_gmdp")) TR_gmdp <- c(475, 975, 1975, 2475, 4975, 9975)
if(!exists("tag")) tag <- "site"  # only one scenario now

# ---------------------------------------------------------------------
# 1–4. Recorrer modelos hermanos de gem_path y procesar cada uno
# ---------------------------------------------------------------------
for (gem_path in list.dirs(dirname(gem_path), full.names = TRUE, recursive = FALSE)) {

  # 1. Localizar carpetas Vs30 y elegir Vref
  stopifnot(dir.exists(gem_path))
  Vs30_gem <- list.dirs(gem_path, full.names = TRUE, recursive = FALSE)
  Vs30_gem <- Vs30_gem[
    sapply(Vs30_gem, function(d) length(list.files(d, all.files = TRUE, no.. = TRUE)) > 0)
  ] |> basename() |> sort()

  if ("760" %in% Vs30_gem) {
    Vref_gmdp <- 760
  } else if ("3000" %in% Vs30_gem) {
    Vref_gmdp <- 3000
  } else {
    stop("Missing Rock Site Conditions model (760 or 3000)")
  }

  # 2. Construir UHS & AEP para este modelo, etiquetando roca con el nombre del modelo
  if (length(Vs30_gem) > 0) {
    LIST <- future_lapply(
      Vs30_gem,
      function(vs) {
        newmark::buildGMDP(
          path = file.path(gem_path, vs),
          vref = vs,
          IDo  = basename(gem_path),   # etiqueta de roca = nombre del modelo
          TRo  = TR_gmdp
        )
      },
      future.seed = TRUE,
      future.scheduling = 4L
    )

    # Tidy por chunk antes de anexar: elimina columnas no usadas y define SaF/AF
    LIST <- lapply(LIST, function(x) {
      x$UHSTable[, `:=`(lat = NULL, lon = NULL, depth = NULL, ITo = NULL, SaF = Sa, AF = 1)]
      x
    })

    UHSTable <- data.table::rbindlist(list(
      UHSTable,
      data.table::rbindlist(lapply(LIST, `[[`, "UHSTable"), use.names = TRUE) |> unique()
    ), use.names = TRUE)

    AEPTable <- data.table::rbindlist(list(
      AEPTable,
      data.table::rbindlist(lapply(LIST, `[[`, "AEPTable"), use.names = TRUE) |> unique()
    ), use.names = TRUE)
  }

  # 3. Subconjunto base en roca para este modelo (criterio físico Vs30==Vref)
  UHSRock <- UHSTable[
    ID == basename(gem_path) &
    TR %in% TR_gmdp &
    Vs30 == Vref_gmdp & Vref == Vref_gmdp,
    .(ID, TR, Tn, Sa, p, Vs30, Vref)
  ]

  # 4. Respuesta de sitio para Vs30 objetivo; etiqueta de sitio = "modelo.tag"
  Vs30_targets <- sort(unique(as.numeric(c(Vs30_gem, Vs30_gmdp))))
  SaFTable <- data.table::rbindlist(
    lapply(Vs30_targets, function(vs) {
      UHSRock[, {
        input <- .SD[, .(Tn, p, Sa)]
        out <- newmark::fitSaF(
          input,
          vs30 = vs,
          vref = Vref_gmdp,
          ns   = 1000,
          Rrup = Rrup_gmdp,
          p_TARGET = c(0.05,0.10,0.16,0.50,0.84,0.90,0.95)
        )
        out[, `:=`(Vs30 = vs, Vref = Vref_gmdp, ID = paste0(basename(gem_path), ".", tag))]
        out
      }, by = .(TR)]
    }),
    use.names = TRUE
  )

  UHSTable <- data.table::rbindlist(list(UHSTable, SaFTable), use.names = TRUE)
}

# ---------------------------------------------------------------------
# 5. Envelope global único en SaF (excluye previamente el propio ID_max)
# ---------------------------------------------------------------------
AUX <- UHSTable[ID != ID_max,
  .SD[which.max(SaF)],
  by = .(TR, Vs30, p, Tn)
]
AUX[, ID := ID_max]
UHSTable <- data.table::rbindlist(list(UHSTable, AUX), use.names = TRUE)

future::plan(sequential)
