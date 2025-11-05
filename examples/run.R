# run.R
# Detect package root (folder containing DESCRIPTION)
.find_pkg_root <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Package root not found from ", start)
    path <- parent
  }
}

pkg_root <- tryCatch(.find_pkg_root(), error = function(e) getwd())

devtools::load_all(pkg_root, quiet = TRUE)

# Paths to examples folder and data
root <- file.path(pkg_root, "examples")

# ======================================================================
# Load example parameters
data_path <- file.path(root, "data")
stopifnot(file.exists(file.path(data_path, "DATA.R")))
source(file.path(data_path, "DATA.R"), local = TRUE)
# ======================================================================

# Build or load UHS/AEP
gem_path  <- file.path(root, "oq", "output", "gem")
  source(file.path(root, "runUHS.R"), local = TRUE)
  saveRDS(UHSTable, file.path(data_path, "UHSTable.Rds"))
  saveRDS(AEPTable, file.path(data_path, "AEPTable.Rds"))

# Build or load Disaggregation (if script exists)
if (!file.exists(file.path(data_path, "DEQTable.Rds"))) {
  if (file.exists(file.path(root, "runDEQ.R"))) {
    source(file.path(root, "runDEQ.R"), local = TRUE)
    saveRDS(DEQTable, file.path(data_path, "DEQTable.Rds"))
  }
} else {
  DEQTable <- readRDS(file.path(data_path, "DEQTable.Rds"))
}

# Build or load ASCE table
if (!file.exists(file.path(data_path, "ASCETable.Rds"))) {
  if (file.exists(file.path(root, "runASCE.R"))) {
    source(file.path(root, "runASCE.R"), local = TRUE)
    saveRDS(ASCETable, file.path(data_path, "ASCETable.Rds"))
  }
} else {
  ASCETable <- readRDS(file.path(data_path, "ASCETable.Rds"))
}

# Build or load shear period table
if (!file.exists(file.path(data_path, "ShearTable.Rds"))) {
  if (file.exists(file.path(root, "runTs.R"))) {
    source(file.path(root, "runTs.R"), local = TRUE)
    saveRDS(ShearTable[], file.path(data_path, "ShearTable.Rds"))
  }
} else {
  ShearTable <- readRDS(file.path(data_path, "ShearTable.Rds"))
}

# Build or load displacement table
if (!file.exists(file.path(data_path, "DnTable.Rds"))) {
  if (file.exists(file.path(root, "runDn.R"))) {
    source(file.path(root, "runDn.R"), local = TRUE)
    saveRDS(DnTable[], file.path(data_path, "DnTable.Rds"))
  }
} else {
  DnTable <- readRDS(file.path(data_path, "DnTable.Rds"))
}

# Build or load kmax table
if (!file.exists(file.path(data_path, "kmaxTable.Rds"))) {
  if (file.exists(file.path(root, "runKmax.R"))) {
    source(file.path(root, "runKmax.R"), local = TRUE)
    saveRDS(kmaxTable, file.path(data_path, "kmaxTable.Rds"))
    # saveRDS(khTable, file.path(data_path, "khTable.Rds"))
  }
} else {
  kmaxTable <- readRDS(file.path(data_path, "kmaxTable.Rds"))
  # khTable <- readRDS(file.path(data_path, "khTable.Rds"))
}


# ======================================================================

Tn_gmdp <- UHSTable[Tn > 0]$Tn |>
  unique() |>
  sort()
p_gmdp <- UHSTable$p |> unique()
Tn_PGA <- min(Tn_gmdp)
TR_gmdp <- UHSTable$TR |>
  unique() |>
  as.numeric() |>
  sort()
Vs30_gmdp <- UHSTable$Vs30 |>
  unique() |>
  as.numeric() |>
  sort()
Vref_gmdp <- UHSTable$Vref[1] |> as.numeric()
# PATCH
#
IDg_gmdp <- DnTable$IDg |> unique()
IDm_gmdp <- DnTable$IDm |> unique()
# Sort sections
IDg_gmdp <- IDg_gmdp[order(as.numeric(sub("^S", "", IDg_gmdp)))]
# Copy files

#file.copy(file.path("../mapper/maps/index.html"),to=file.path(file.path(root,"map/index.html")),overwrite = TRUE)
#file.copy(file.path("../mapper/maps/epicenters.csv"),to=file.path(file.path(root,"map/epicenters.csv")),overwrite = TRUE)