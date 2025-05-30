### DESCRIPTION ▸ `Description:` field

```
Implements a unified Monte-Carlo framework for propagating (i) the epistemic
uncertainty of uniform-hazard spectral ordinates Sa(Tn)—captured through their
quantiles extracted from hazard curves—and (ii) the aleatory variability of
empirical ground-failure models expressed as log-normal random variables.
Core builders:

  • fitDn()  – combines multiple sliding-block Newmark-displacement models
               (AM88, YG91, JB07, BT07, SR08, BM17, BM19) with sampled Sa, PGA,
               PGV, AI, period-dependent scaling and user-defined model weights
               to produce scenario-specific displacement quantiles.

  • fitSaF() – convolves Sa(Tn) with the Stewart & Seyhan (2017) site-factor
               FST17(Tn,PGA,Vs30) to obtain amplified spectral ordinates for
               site-response analysis.

  • sampleSa() – draws Monte-Carlo realisations of Sa at arbitrary periods
                 consistent with user-supplied Sa(Tn,q) quantiles via the sdQ
                 method.

The package thus delivers end-to-end, reproducible uncertainty convolution for
probabilistic slope-stability and site-response studies, using only
data.table and base-R dependencies.
```

---

### README.md  (minimal working example)

````markdown
# newmark <img src="man/figures/logo.png" align="right" height="140"/>

Monte-Carlo convolution of **uniform-hazard spectra** with
**empirical ground-failure models** for:

* **Sliding-block Newmark displacements** (multiple models, weighted)
* **Non-linear site factors** (Stewart & Seyhan 2017)

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("your-github/newmark")
````

---

## Quick Start

```r
library(newmark)
library(data.table)

## ── 1  Build a toy UHS (Tn, Sa, p = c(0.05, 0.5, 0.95, "mean")) ──────────
UHS <- data.table(
  Tn = rep(c(0.01, 0.1, 0.2, 0.5, 1, 2), each = 4),
  p  = rep(c(0.05, 0.5, 0.95, "mean"), 6),
  Sa = c(0.25, 0.35, 0.55, 0.38,      # 0.01 s ≈ PGA
         0.38, 0.55, 0.88, 0.60,      # 0.10 s
         0.42, 0.65, 1.10, 0.72,      # …
         0.30, 0.48, 0.82, 0.55,
         0.19, 0.31, 0.55, 0.36,
         0.08, 0.14, 0.26, 0.17)
)

## ── 2  Newmark-displacement quantiles with default model weights ─────────
fitDn(
  UHS    = UHS,
  ky     = 0.15,    # yield acceleration (g)
  Ts     = 1.0,     # period of sliding mass (s)
  Mw     = 6.5,     # scenario magnitude
  NS     = 1000     # Monte-Carlo samples per model
)

#       p        Dn
#    <num>     <cm>
# 1: 0.05  0.21
# 2: 0.50  2.47
# 3: 0.95 14.62
# 4: mean  3.35

## ── 3  Site-amplified SaF(Tn) quantiles at Vs30 = 400 m/s ───────────────
fitSaF(
  UHS   = UHS,
  vs30  = 400,
  Tn    = 0.5,    # period of interest
  NS    = 1000
)

#       p    SaF
# 1: 0.05 0.47
# 2: 0.50 0.69
# 3: 0.95 1.05
# 4: mean 0.72
```

---

## Function Reference (most common)

| Function     | Purpose                                                      | Key inputs                    |
| ------------ | ------------------------------------------------------------ | ----------------------------- |
| `sampleSa()` | Draw Sa realisations at any period based on UHS quantiles    | `UHS`, `Td`, `n`              |
| `fitDn()`    | Weighted Newmark-displacement quantiles from multiple models | `UHS`, `ky`, `Ts`, `Mw`, `NS` |
| `fitSaF()`   | Sa × F<sub>ST17</sub> quantiles for a single Vs30            | `UHS`, `vs30`, `Tn`, `NS`     |

All helpers return **data.table** objects for seamless downstream analysis.

---

## Citation

If you use **newmark** in academic work, please cite the primary reference
of each empirical model you rely on (Ambraseys & Menu 1988, Yegian et al. 1991,
Jibson 2007, Bray & Travasarou 2007, Saygili & Rathje 2008,
Bray & Macedo 2017/2019, Stewart & Seyhan 2017) and this package.

---

*© 2025 Alejandro Verri Kozlowski*

