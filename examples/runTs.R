# nolint start
root <- here::here()
if(!exists("NTS")){NTS <- 10}
if(!exists("Hs")){
    Hs <- seq(Hs.min, Hs.max, length.out = NTS)
}


lo <- seq(lo.min, lo.max, length.out = length(Hs))
s <- seq(s.min, s.max, length.out = length(Hs))
b <- (2 * Hs * s / (1 / lo - 1))
bmax <- (b + 2 * Hs * s)
IDg <- paste0("S", round(Hs, 0))
GeometryTable <- data.table::data.table(
    IDg = IDg,
    b = b |> round(0),
    bmax = bmax |> round(0),
    s = s |> round(2),
    Hs = Hs |> round(0),
    beta = atan(1 / s) * 180 / pi
)
p_TARGET <- c(0.02, 0.05, 0.10, 0.16, "mean", 0.84, 0.90, 0.95, 0.98)

MaterialTable <- GeometryTable[CJ(IDm, IDg = GeometryTable$IDg), on = .(IDg), allow.cartesian = TRUE]
MaterialTable[, Vref := Vref_gmdp[1]] # Add Reference bedrock

# Calculate
ShearTable <- MaterialTable[,
    {
        OUT <- dsra::getSiteProperties(Hs = Hs, USCS = uscs[[IDm]], h = 1.00, NR = 100, levels = p_TARGET, Vref = Vref)
        OUT[, Hs := NULL] # remove the duplicate
    },
    by = .(Hs, b, s, Vref, IDm, IDg)
][]



# Characteristic Roots
ShearTable[, lo := round((b / (2 * Hs * s + b)), 3)]
ShearTable[, an := mapply(function(mo, lo) {
    dsra::getCylinderRoots(mo = mo, lo = lo, no = 1, model = "nlm", extrapolate = TRUE) |>
        round(4) |>
        suppressWarnings()
}, mo, lo)]
# Fundamental periods ----
ShearTable[, Ts := (4 * pi * Hs / (an * (2 - mo) * VSo)) |> round(3)]
# _slides/ is the home folder. up one level to data/
