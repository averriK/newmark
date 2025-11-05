# runDEQ.R
# nolint start
root <- here::here()

# NOTE: Only ANCOLD prescribes 84 % for Extreme categories.
# Since we don’t have explicit MCE results, we adopt 84 % for all Standards
# for Extreme categories.

# ---------------- 1. ANCOLD --------------------------------------------------
ANCOLD <- data.table(
    Standard = "ANCOLD",
    Class = c("Low", "Significant", "High C", "High B", "High A", "Extreme"),
    Operation = rep(975, 6), # actually is mid‑range OBE (1/475‑1/1 000)
    Closure = c(975, 975, 1975, 4975, 9975, 9975),
    `Post-Closure` = c(975, 975, 1975, 4975, 9975, 9975)
)

# ---------------- 2. CDA (patched) -------------------------------------------
#  • Low Closure & Post‑Closure: use 975
#  • Significant Closure & Post‑Closure: use 2 475
#  • High / Very High / Extreme remain at 9 975 (≈10 000) for passive‑care
CDA <- data.table(
    Standard = "CDA",
    Class = c("Low", "Significant", "High", "Very High", "Extreme"),
    Operation = c(975, 975, 2475, 4975, 9975),
    Closure = c(975, 2475, 9975, 9975, 9975),
    `Post-Closure` = c(975, 2475, 9975, 9975, 9975)
)

# ---------------- 3. GISTM (Low kept at 475) ---------------------------------
GISTM <- data.table(
    Standard = "GISTM",
    Class = c("Low", "Significant", "High", "Very High", "Extreme"),
    Operation = c(475, 975, 2475, 4975, 9975), # 1/200 mapped to 475 AEP slot
    Closure = c(475, 975, 2475, 4975, 9975),
    `Post-Closure` = rep(9975, 5)
)

longify <- function(dt) {
    # 'melt' everything except Standard, Class
    melt(
        dt,
        id.vars       = c("Standard", "Class"),
        variable.name = "Stage",
        value.name    = "TR"
    )
}

# ---------- 4. Combine all tables ---------------------------------------
STDT <- rbindlist(
    list(
        longify(ANCOLD),
        longify(CDA),
        longify(GISTM)
    ),
    use.names = TRUE
)

Standard_TARGET <- unique(STDT$Standard)
Stage_TARGET <- unique(STDT$Stage)
p_TARGET <- "mean"

DEQTable <- data.table()
for (S in Standard_TARGET) {
    for (L in Stage_TARGET) {
        AUX <- UHSTable[STDT[Standard == S & Stage == L], on = .(TR)][p %in% p_TARGET]
        DEQTable <- list(DEQTable, AUX) |> rbindlist()
    }
}
DEQTable[, Stage := as.character(Stage)]
