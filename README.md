# newmark

**Monte-Carlo Ensemble Newmark Displacements**

[![R Version](https://img.shields.io/badge/R-%3E%3D%204.2.0-blue)](https://www.r-project.org/) [![Version](https://img.shields.io/badge/version-0.4.0-green)](https://github.com/averriK/newmark)

R package for probabilistic seismic displacement analysis using Newmark sliding-block method with epistemic uncertainty propagation.

## What is it?

newmark propagates uncertainty from uniform-hazard spectral accelerations through empirical sliding-block and non-linear site-factor models using weighted Monte-Carlo sampling. Combines multiple displacement prediction models with epistemic weighting to deliver scenario-specific displacement quantiles.

## Features

- **PSHA import**: OpenQuake and user-defined hazard curves
- **Site amplification**: Stewart & Seyhan (2017) with spectral correlation
- **Ensemble displacement models**: AM88, YG91, JB07, SR08, BT07, BM17, BM19
- **Spectral correlation**: Baker & Jayaram (2009) for Sa(T) sampling
- **Parallel processing**: Future-based workflows for large parameter grids
- **Epistemic weighting**: Model uncertainty quantification

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

## Documentation

See function documentation via R help:

```r
?newmark
?buildGMDP
?fitSaF
?fitDn
```

Full API: `buildGMDP()`, `fitSaF()`, `fitDn()`, displacement models (7 models: AM88, YG91, JB07, SR08, BT07, BM17, BM19), `rhoBJ()`, `sampleSaCorr()`, `kmaxMC()`, `checkUHS()`, `designUHS()`, `Vs30toSID()`, `approx.spline()`

## Dependencies

- R (>= 4.2)
- data.table, Hmisc, mvtnorm, stats
- future, future.apply (parallel processing)
- readxl, stringr, devtools

## References

Baker, J. W., & Jayaram, N. (2008). Correlation of spectral acceleration values from NGA ground motion models. *Earthquake Spectra*, 24(1), 299-317.

Stewart, J. P., & Seyhan, E. (2013). Semi-empirical nonlinear site amplification and its application in NEHRP site factors. PEER Report 2013/13.

Bray, J. D., & Macedo, J. (2019). Procedure for estimating shear-induced seismic slope displacement for shallow crustal earthquakes. *Journal of Geotechnical and Geoenvironmental Engineering*, 145(12), 04019106.

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
