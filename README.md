# newmark

**Monte-Carlo Ensemble Newmark Displacements**

> **Last updated:** November 13, 2025

R package for probabilistic seismic displacement analysis using Newmark sliding-block method with epistemic uncertainty propagation.

[![R Version](https://img.shields.io/badge/R-%3E%3D%204.2.0-blue)](https://www.r-project.org/)
[![Version](https://img.shields.io/badge/version-0.4.0-green)](https://github.com/averriK/newmark)

## What is it?

newmark propagates uncertainty from uniform-hazard spectral accelerations through empirical sliding-block and non-linear site-factor models using weighted Monte-Carlo sampling. Combines multiple displacement prediction models with epistemic weighting to deliver scenario-specific displacement quantiles.

## Features

- **PSHA import**: OpenQuake and user-defined hazard curves
- **Site amplification**: Stewart & Seyhan (2017) with spectral correlation
- **Ensemble displacement models**: AM88, YG91, JB07, SR08, BT07, BM17, BM19
- **Spectral correlation**: Baker & Jayaram (2009) for Sa(T) sampling
- **Parallel processing**: Future-based workflows for large parameter grids

## Installation

```r
devtools::install_github("averriK/newmark")
```

## Usage

```r
library(newmark)

# Import PSHA
GMDP <- buildGMDP(
  path = "openquake/output/760",
  IDo = "gem",
  vref = 760,
  TRo = c(475, 2475, 4975)
)

# Site amplification
SaFTable <- fitSaF(
  uhs = GMDP$UHSTable,
  vs30 = 270,
  vref = 760,
  ns = 1000,
  Rrup = 100
)

# Displacement analysis
DnTable <- fitDn(
  uhs = SaFTable,
  ky = 0.15,
  Ts = 0.25,
  Mw = 6.5,
  NS = 200,
  Rrup = 100
)
```

## API overview

- buildGMDP(path, IDo="GEM", engine="openquake"|"user", vref=760, TRo=...): Import hazard (zip from OpenQuake or user AEP.xlsx) and build AEP/UHS tables (plus disaggregation if available).
- fitSaF(uhs, vs30, vref=760, ns=1000, Rrup=100): Site amplification (Stewart & Seyhan 2017) with spectral correlation; returns SaF and AF by Tn and p.
- fitDn(uhs, ky, Ts, Mw=6.5, NS=30, Rrup=100): Ensemble Newmark displacements with epistemic weights; returns quantiles and mean.
- Displacement models: Dn_AM88, Dn_YG91, Dn_JB07, Dn_SR08, Dn_BT07, Dn_BM17, Dn_BM19 (log-mean/log-sd on cm or ln scale per model).
- Helpers: buildQSpline, interpolateSaTable, rhoBJ, sampleSaCorr, kmaxMC, checkUHS, designUHS.
- Utilities: Vs30toSID, approx.spline.

Note: For engine="user" in buildGMDP(), place an AEP.xlsx file under the provided path with the expected columns.

## Examples

Complete workflows available in:

- `~/github/psha/R/runDn.R` - Displacement analysis workflow
- `~/github/psha/R/runKmax.R` - Seismic coefficient analysis

## Documentation

Function documentation available via R help:

```r
?newmark
?buildGMDP
?fitSaF
?fitDn
```

## Dependencies

- R (>= 4.2)
- data.table, Hmisc, mvtnorm
- future, future.apply
- readxl, stringr, devtools

## License

Custom license - see [LICENSE](LICENSE)

## Citation

```bibtex
@software{newmark2024,
  author = {Verri Kozlowski, Alejandro},
  title = {newmark: Monte-Carlo Ensemble Newmark Displacements},
  year = {2024},
  version = {0.4.0},
  url = {https://github.com/averriK/newmark}
}
```

---

**Author:** Alejandro Verri Kozlowski  
**Email:** averri@fi.uba.ar  
**ORCID:** [0000-0002-8535-1170](https://orcid.org/0000-0002-8535-1170)
