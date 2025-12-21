# newmark

**Monte-Carlo Ensemble Newmark Displacements**

[![R Version](https://img.shields.io/badge/R-%3E%3D%204.2.0-blue)](https://www.r-project.org/) [![Version](https://img.shields.io/badge/version-0.4.0-green)](https://github.com/averriK/newmark)

Propagates uncertainty from uniform-hazard spectral accelerations through empirical sliding-block and non-linear site-factor models using weighted Monte-Carlo sampling. Combines multiple displacement prediction models with epistemic weighting to deliver scenario-specific displacement quantiles.

## Contents

- [Overview](#overview)
- [Data conventions](#data-conventions)
- [Core workflows](#core-workflows)
- [Installation](#installation)
- [Exported API](#exported-api)
- [Dependencies](#dependencies)
- [References](#references)
- [License](#license)
- [Author](#author)

## Overview

The package implements a three-stage workflow for probabilistic seismic displacement analysis:

1. **Hazard import** (`buildGMDP`): assembles annual-exceedance-probability (AEP) tables, uniform-hazard spectra (UHS), and optional magnitude-distance disaggregation from OpenQuake or user-supplied files.

2. **Site amplification** (`fitSaF`): applies the Stewart & Seyhan (2017) non-linear site factor model to rock spectra, with spectral correlation via Baker & Jayaram (2009) when quantile data is available.

3. **Displacement estimation** (`fitDn`): samples correlated spectral ordinates and runs an epistemic ensemble of seven Newmark displacement models (AM88, YG91, JB07, SR08, BT07, BM17, BM19) to produce weighted displacement quantiles.

All outputs are returned as `data.table` objects.

## Data conventions

- **Periods**: `Tn = 0` denotes PGA; `Tn = -1` denotes PGV (if present).
- **Probabilities**: quantile levels are numeric (e.g., `0.16`, `0.50`, `0.84`); the mean is labelled `"mean"`.
- **Units**: spectral acceleration in g, displacement in cm, velocity in cm/s, periods in s, Vs30 in m/s.
- **Hazard tables**: columns include `Tn`, `p`, `Sa`, `TR` (return period in years), `AEP` (annual exceedance probability), `POE` (probability of exceedance).
- **UHS tables**: uniform-hazard spectra interpolated to a regular TR grid; columns include `Tn`, `p`, `Sa`, `TR`, `Vs30`, `Vref`, `SaF` (site-amplified Sa), `AF` (amplification factor).

## Core workflows

### 1. Import PSHA hazard data

```r
GMDP <- buildGMDP(
  path = "openquake/output/760",
  IDo = "gem",
  engine = "openquake",
  vref = 760,
  TRo = c(475, 2475, 4975)
)
```

**Inputs**: directory containing OpenQuake `.zip` files or user-supplied `AEP.xlsx`.

**Outputs**: list with `AEPTable`, `UHSTable`, and `RMwTable` (disaggregation, if available).

**Function**: `buildGMDP()` (R/buildGMDP.R:20–62)

### 2. Site amplification with spectral correlation

```r
SaFTable <- fitSaF(
  uhs = GMDP$UHSTable,
  vs30 = 270,
  vref = 760,
  ns = 1000,
  Rrup = 100
)
```

**Inputs**: uniform-hazard spectrum with columns `Tn`, `p`, `Sa`; target `vs30` and reference `vref`.

**Outputs**: `data.table` with columns `Tn`, `p`, `Sa`, `SaF`, `AF`.

**Modes**:
- When input contains only `p == "mean"`, rock Sa is treated as deterministic and dispersion comes solely from the site factor model.
- When input contains numeric quantiles, ln(Sa) is sampled via Gaussian copula with Baker & Jayaram (2009) correlation between PGA and Sa(Tn).

**Function**: `fitSaF()` (R/fitSaF.R:35–155)

### 3. Ensemble displacement analysis

```r
DnTable <- fitDn(
  uhs = SaFTable,
  ky = 0.15,
  Ts = 0.25,
  Mw = 6.5,
  NS = 200,
  Rrup = 100
)
```

**Inputs**: uniform-hazard spectrum with columns `Tn`, `p`, `Sa`; yield acceleration `ky` (g), sliding-mass period `Ts` (s), moment magnitude `Mw`, Monte-Carlo sample size `NS`, rupture distance `Rrup` (km).

**Outputs**: `data.table` with columns `p`, `Dn` (displacement in cm).

**Process**: samples correlated spectral ordinates at four periods (PGA, 1.0·Ts, 1.3·Ts, 1.5·Ts) using a Gaussian copula; runs seven displacement models with epistemic weights (YG91=1/19, AM88=2/19, JB07=2/19, BT07=3/19, SR08=3/19, BM17=4/19, BM19=4/19); returns weighted quantiles.

**Function**: `fitDn()` (R/fitDn.R:16–137)

### 4. Invert displacement to maximum yield acceleration

```r
kmaxTable <- kmaxMC(
  Dn = DnTable$Dn,
  ky = DnTable$ky,
  p = DnTable$p,
  Da = c(2.5, 25, 250),
  p_TARGET = c(0.16, "mean", 0.84),
  ns = 5e4
)
```

**Inputs**: displacement values `Dn`, yield accelerations `ky`, probability labels `p`, target displacements `Da` (cm), output probability levels `p_TARGET`, Monte-Carlo sample size `ns`.

**Outputs**: `data.table` with columns `Da`, `p`, `kmax`.

**Function**: `kmaxMC()` (R/kmaxMC.R:21–109)

## Installation

```r
# Install from GitHub
devtools::install_github("averriK/newmark")
```

## Exported API

The package exports 23 functions (NAMESPACE:1–29):

### Main workflows
- `buildGMDP()` — import PSHA hazard data (OpenQuake or user files)
- `fitSaF()` — site amplification with spectral correlation (Stewart & Seyhan 2017)
- `fitDn()` — Monte-Carlo ensemble displacement analysis

### Displacement models (7 models)
- `Dn_AM88()` — Ambraseys & Menu (1988)
- `Dn_YG91()` — Yegian et al. (1991)
- `Dn_JB07()` — Jibson (2007)
- `Dn_SR08()` — Saygili & Rathje (2008)
- `Dn_BT07()` — Bray & Travasarou (2007)
- `Dn_BM17()` — Bray & Macedo (2017)
- `Dn_BM19()` — Bray & Macedo (2019)

### Site factor models
- `F_ST17()` — Stewart & Seyhan (2017) non-linear site factor
- `SaF_ST17()` — canonical wrapper for site-amplified Sa

### Spectral correlation and sampling
- `rhoBJ()` — Baker & Jayaram (2009) correlation model
- `sampleSaCorr()` — correlated Sa sampling via Gaussian copula
- `buildQSpline()` — monotone quantile spline for ln(Sa)

### Utilities
- `kmaxMC()` — invert Dn(ky) to obtain k_max distribution
- `interpolateSaTable()` — log-log interpolation of Sa at arbitrary period
- `checkUHS()` — validate UHS monotonicity and detect duplicates
- `designUHS()` — ASCE 7-22 design spectrum (MCER or 2/3·MCER)
- `Vs30toSID()` — map Vs30 to site class identifier
- `approx.spline()` — spline interpolation wrapper

## Dependencies

R (>= 4.2)

**Imports**: data.table, Hmisc, mvtnorm, stats, readxl, stringr, future, future.apply, devtools

**Suggests**: testthat (>= 3.0.0), knitr, rmarkdown

(DESCRIPTION:21–36)

## References

Ambraseys, N. N., & Menu, J. M. (1988). Earthquake-induced ground displacements. *Earthquake Engineering & Structural Dynamics*, 16(7), 985–1006.

Baker, J. W., & Jayaram, N. (2009). Correlation model for spatially distributed ground-motion intensities. *Earthquake Engineering & Structural Dynamics*, 38(8), 951–972.

Bray, J. D., & Macedo, J. (2017). Procedure for estimating shear-induced seismic slope displacement for subduction zone earthquakes. *Journal of Geotechnical and Geoenvironmental Engineering*, 143(12), 04017106.

Bray, J. D., & Macedo, J. (2019). Procedure for estimating shear-induced seismic slope displacement for shallow crustal earthquakes. *Journal of Geotechnical and Geoenvironmental Engineering*, 145(12), 04019106.

Bray, J. D., & Travasarou, T. (2007). Simplified procedure for estimating earthquake-induced deviatoric slope displacements. *Journal of Geotechnical and Geoenvironmental Engineering*, 133(4), 381–392.

Jibson, R. W. (2007). Regression models for estimating coseismic landslide displacement. *Engineering Geology*, 91(2-4), 209–218.

Saygili, G., & Rathje, E. M. (2008). Empirical predictive models for earthquake-induced sliding displacements of slopes. *Journal of Geotechnical and Geoenvironmental Engineering*, 134(6), 790–803.

Stewart, J. P., & Seyhan, E. (2017). Semi-empirical nonlinear site amplification and its application in NEHRP site factors. PEER Report 2013/13.

Yegian, M. K., Marciano, E. A., & Ghahraman, V. G. (1991). Earthquake-induced permanent deformations: probabilistic approach. *Journal of Geotechnical Engineering*, 117(1), 35–50.

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

## Author

**Alejandro Verri Kozlowski**  
Email: averri@fi.uba.ar  
ORCID: [0000-0002-8535-1170](https://orcid.org/0000-0002-8535-1170)  
Affiliation: Universidad de Buenos Aires, Facultad de Ingeniería
