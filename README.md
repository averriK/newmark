# newmark

**Monte-Carlo Ensemble Newmark Displacements**

R package for probabilistic seismic displacement analysis of earth structures using the Newmark sliding-block method. Propagates uncertainty from uniform-hazard spectral accelerations through empirical sliding-block and non-linear site-factor models using weighted Monte-Carlo sampling.

[![R Version](https://img.shields.io/badge/R-%3E%3D%204.2.0-blue)](https://www.r-project.org/)
[![Version](https://img.shields.io/badge/version-0.4.0-green)](https://github.com/averriK/newmark)

## Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Core Functions](#core-functions)
- [Workflow](#workflow)
- [Theory and Mathematical Background](#theory-and-mathematical-background)
  - [Site Amplification: Stewart & Seyhan (2017)](#site-amplification-stewart--seyhan-2017)
  - [Newmark Displacement Models](#newmark-displacement-models)
  - [Epistemic Weighting](#epistemic-weighting)
  - [Monte Carlo Sampling and Spectral Correlation](#monte-carlo-sampling-and-spectral-correlation)
- [Dependencies](#dependencies)
- [References](#references)

---

## Installation

```r
# Install from GitHub
# install.packages("devtools")
devtools::install_github("averriK/newmark")
```

## Quick Start

### 1. Build uniform hazard spectra from OpenQuake output

```r
library(newmark)
library(data.table)

# Load PSHA hazard curves (OpenQuake zip files)
GMDP <- buildGMDP(
  path = "path/to/openquake/output/760",
  IDo = "gem",
  engine = "openquake",
  vref = 760,  # Reference Vs30 (m/s)
  TRo = c(475, 975, 1975, 2475, 4975, 9975)
)

UHSTable <- GMDP$UHSTable  # Uniform hazard spectra
AEPTable <- GMDP$AEPTable  # Annual exceedance probabilities
```

### 2. Compute site amplification for target Vs30

```r
# Extract rock spectrum for a specific return period
UHSRock <- UHSTable[
  ID == "gem" & TR == 2475 & Vs30 == 760 & Vref == 760,
  .(Tn, p, Sa)
]

# Amplify to soft soil site (Vs30 = 270 m/s)
SaFTable <- fitSaF(
  uhs = UHSRock,
  vs30 = 270,
  vref = 760,
  ns = 1000,    # Monte Carlo sample size
  Rrup = 100,   # Rupture distance (km) for spectral correlation
  p_TARGET = c(0.05, 0.10, 0.16, 0.50, 0.84, 0.90, 0.95)
)

# SaFTable contains: Tn, p, Sa (rock), SaF (site-amplified), AF (amplification factor)
```

### 3. Compute Newmark displacements

```r
# Compute displacement for specified yield acceleration and period
DnTable <- fitDn(
  uhs = SaFTable,
  ky = 0.15,     # Yield acceleration (g)
  Ts = 0.25,     # Fundamental period of sliding mass (s)
  Mw = 6.5,      # Moment magnitude
  NS = 200,      # Monte Carlo sample size
  Rrup = 100     # Rupture distance (km)
)

# DnTable contains: p (quantile), Dn (displacement in cm)
print(DnTable)
```

---

## Core Functions

### PSHA Import and Processing

| Function | Description | Input | Output |
|----------|-------------|-------|--------|
| `buildGMDP` | Import OpenQuake or user hazard curves | Hazard curve files | UHSTable, AEPTable, RMwTable |

**Key parameters:**
- `path`: Directory with hazard data files
- `engine`: "openquake" or "user"
- `vref`: Reference Vs30 (m/s), typically 760 or 3000
- `TRo`: Target return periods (years)

### Site Response

| Function | Description | Input | Output |
|----------|-------------|-------|--------|
| `fitSaF` | Site amplification (Stewart & Seyhan 2017) | Rock UHS, target Vs30 | Site-amplified spectra |
| `F_ST17` | Evaluate site factor F_ST17 | PGA, Tn, vs30, vref | muLnF, sdLnF |

**fitSaF modes:**
- Mode A: Input contains only p = "mean" → dispersion from F_ST17 aleatory uncertainty
- Mode B: Input contains quantiles (0.05, 0.16, ...) → sample rock Sa with Baker & Jayaram correlation

### Newmark Displacements

| Function | Description | Input | Output |
|----------|-------------|-------|--------|
| `fitDn` | Ensemble displacement model | Site spectrum, ky, Ts, Mw | Displacement quantiles |

**Individual displacement models:**
- `Dn_AM88`: Ambraseys & Menu (1988) rigid block
- `Dn_YG91`: Yegian et al. (1991) rigid block
- `Dn_JB07`: Jibson (2007) with Arias intensity
- `Dn_SR08`: Saygili & Rathje (2008) with Arias intensity
- `Dn_BT07`: Bray & Travasarou (2007) flexible block
- `Dn_BM17`: Bray & Macedo (2017) subduction events
- `Dn_BM19`: Bray & Macedo (2019) shallow-crustal with pulse effects

### Utilities

| Function | Description |
|----------|-------------|
| `sampleSaCorr` | Sample correlated spectral ordinates |
| `rhoBJ` | Baker & Jayaram (2009) correlation coefficient |
| `buildQSpline` | Empirical quantile function (Hyman spline) |
| `interpolateSaTable` | Log-log interpolation for missing periods |

---

## Workflow

The typical analysis chain:

```r
library(newmark)
library(future)

# Enable parallel processing
plan(multisession, workers = parallel::detectCores() - 1)

# 1. Build UHS from PSHA
GMDP <- buildGMDP(
  path = "examples/oq/output/gem/760",
  IDo = "gem",
  vref = 760,
  TRo = c(475, 975, 1975, 2475, 4975, 9975)
)

UHSTable <- GMDP$UHSTable

# 2. Site amplification for multiple Vs30
Vs30_targets <- c(180, 270, 360, 560, 760)

SaFTable <- rbindlist(lapply(Vs30_targets, function(vs) {
  UHSRock <- UHSTable[Vs30 == 760 & Vref == 760, .(TR, Tn, p, Sa)]

  UHSRock[, fitSaF(
    .SD[, .(Tn, p, Sa)],
    vs30 = vs,
    vref = 760,
    ns = 1000,
    Rrup = 100
  ), by = .(TR)]
}))

# 3. Displacement analysis for parameter grid
ky_vec <- seq(0.05, 0.50, length.out = 20)
Ts_vec <- c(0.15, 0.25, 0.35, 0.50)
TR_vec <- c(475, 975, 2475, 4975, 9975)

DnTable <- rbindlist(lapply(TR_vec, function(tr) {
  rbindlist(lapply(Ts_vec, function(ts) {
    rbindlist(lapply(ky_vec, function(k) {
      UHS <- SaFTable[TR == tr & Vs30 == 270, .(Tn, p, Sa = SaF)]

      out <- fitDn(
        uhs = UHS,
        ky = k,
        Ts = ts,
        Mw = 6.5,
        NS = 200,
        Rrup = 100
      )

      out[, `:=`(TR = tr, ky = k, Ts = ts, Vs30 = 270)]
    }))
  }))
}))

# 4. Save results
fwrite(DnTable, "DnTable.csv")

# Restore sequential processing
plan(sequential)
```

---

## Theory and Mathematical Background

### Site Amplification: Stewart & Seyhan (2017)

The site amplification factor $F$ relates rock spectral acceleration $\text{Sa}(T_n, V_{\text{ref}})$ to site-amplified $\text{SaF}(T_n, V_{s30})$:

$$
\ln[\text{SaF}(T_n)] = \ln[\text{Sa}(T_n, V_{\text{ref}})] + \ln F(T_n, \text{PGA}, V_{s30})
$$

where the site factor is decomposed into:

$$
\ln F = f_L(T_n, V_{s30}) + f_{NL}(T_n, \text{PGA}, V_{s30}) + f_I(T_n, V_{\text{ref}}) + \varepsilon
$$

#### Linear Component: $f_L(T_n, V_{s30})$

Period-dependent linear amplification:

$$
f_L = \begin{cases}
c(T_n) \cdot \ln\left(\frac{V_1(T_n)}{760}\right) & \text{if } V_{s30} \le V_1(T_n) \\
c(T_n) \cdot \ln\left(\frac{V_{s30}}{760}\right) & \text{if } V_1(T_n) < V_{s30} \le V_2(T_n) \\
c(T_n) \cdot \ln\left(\frac{V_2(T_n)}{760}\right) + \frac{c(T_n)}{2} \cdot \ln\left(\frac{V_{s30}}{V_2(T_n)}\right) & \text{if } V_{s30} > V_2(T_n)
\end{cases}
$$

where $c(T_n)$, $V_1(T_n)$, $V_2(T_n)$ are period-dependent coefficients interpolated from tabulated values in F_ST17.R:124-167.

#### Non-Linear Component: $f_{NL}(T_n, \text{PGA}, V_{s30})$

PGA-dependent reduction for soft sites:

$$
f_{NL} = \begin{cases}
a(T_n, V_{s30}) \cdot \ln\left(\frac{\exp(\ln \text{PGA})}{b(T_n)} + 1\right) & \text{if } V_{s30} < V_c(T_n) \\
0 & \text{otherwise}
\end{cases}
$$

where:

$$
\begin{aligned}
a(T_n, V_{s30}) &= f_4(T_n) \cdot \left[\exp(f_5(T_n) \cdot (\min(V_{s30}, 3000) - 360)) - \exp(f_5(T_n) \cdot (3000 - 360))\right] \\
b(T_n) &= C_{760/3000} \cdot f_3(T_n) \\
C_{760/3000} &= \begin{cases} 2.275 & \text{if } V_{\text{ref}} = 760 \\ 1.0 & \text{if } V_{\text{ref}} = 3000 \end{cases}
\end{aligned}
$$

Functions $f_3$, $f_4$, $f_5$, $V_c$ are interpolated from tables in F_ST17.R:263-303.

#### Impedance Adjustment: $f_I(T_n, V_{\text{ref}})$

Correction when reference velocity is hard rock:

$$
f_I = \begin{cases}
F_{760}(T_n) & \text{if } V_{\text{ref}} = 3000 \\
0 & \text{if } V_{\text{ref}} = 760
\end{cases}
$$

where $F_{760}(T_n)$ is interpolated from F_ST17.R:201-233.

#### Aleatory Uncertainty: $\varepsilon \sim N(0, \sigma_{\ln F})$

Total standard deviation:

$$
\sigma_{\ln F} = \sqrt{\sigma_L^2 + \sigma_{NL}^2 + \sigma_I^2}
$$

**Linear component $\sigma_L$:**

$$
\sigma_L = \begin{cases}
\sigma_{L,I} - 2 \cdot (\sigma_{L,I} - \sigma_{Vc}) \cdot \frac{V_{s30} - V_l}{V_f - V_l} + (\sigma_{L,I} - \sigma_{Vc}) \cdot \left(\frac{V_{s30} - V_l}{V_f - V_l}\right)^2 & \text{if } V_{s30} \le V_f(T_n) \\
\sigma_{Vc}(T_n) & \text{if } V_f(T_n) < V_{s30} \le V_2(T_n) \\
\sigma_{Vc}(T_n) + (\sigma_{U,I} - \sigma_{Vc}) \cdot \left(\frac{V_{s30} - V_2(T_n)}{V_u - V_2(T_n)}\right)^2 & \text{if } V_{s30} > V_2(T_n)
\end{cases}
$$

where $V_l = 200$ m/s, $V_u = 2000$ m/s, and $\sigma_{L,I}$, $\sigma_{Vc}$, $\sigma_{U,I}$, $V_f$ are interpolated from F_ST17.R:168-200.

**Non-linear component $\sigma_{NL}$:**

$$
\sigma_{NL} = \begin{cases}
\sigma_c(T_n) & \text{if } V_{s30} < 300 \\
-\sigma_c(T_n) \cdot \frac{\ln(V_{s30}/300)}{\ln(1000/300)} + \sigma_c(T_n) & \text{if } 300 \le V_{s30} < 1000 \\
0 & \text{if } V_{s30} \ge 1000
\end{cases}
$$

where $\sigma_c(T_n)$ is interpolated from F_ST17.R:297-303.

**Impedance component $\sigma_I$:**

$$
\sigma_I = \begin{cases}
\sigma_{760}(T_n) & \text{if } V_{\text{ref}} = 3000 \\
0 & \text{if } V_{\text{ref}} = 760
\end{cases}
$$

where $\sigma_{760}(T_n)$ is interpolated from F_ST17.R:234-262.

**Implementation:** F_ST17.R:108-356

---

### Newmark Displacement Models

The package implements seven empirical displacement models combined using epistemic weights. All models return $\mu_{\ln D_n}$ and $\sigma_{\ln D_n}$ for lognormal displacement distribution.

#### Epistemic Weighting

fitDn.R:64-68:

```r
weight <- data.table(
  ID = c("YG91", "AM88", "JB07", "BT07", "SR08", "BM17", "BM19"),
  w  = c(1, 2, 2, 3, 3, 4, 4)
)
weight[, w := w / sum(w)]  # Normalized: 1/19, 2/19, 2/19, 3/19, 3/19, 4/19, 4/19
```

#### Rigid Sliding-Block Models

##### Dn_AM88: Ambraseys & Menu (1988)

Dn_models.R:8-19:

$$
\begin{aligned}
r &= \min\left(\frac{k_y}{\text{PGA}}, 0.9999\right) \\
\log_{10} D_n &= 0.90 + \log_{10}\left[(1 - r)^{2.53} \cdot r^{-1.09}\right] \\
\sigma_{\log_{10} D_n} &= 0.30
\end{aligned}
$$

Conversion to natural log: $\mu_{\ln D_n} = \mu_{\log_{10} D_n} \cdot \ln(10)$, $\sigma_{\ln D_n} = \sigma_{\log_{10} D_n} \cdot \ln(10)$.

##### Dn_YG91: Yegian et al. (1991)

Dn_models.R:31-46:

$$
\begin{aligned}
r &= \frac{k_y}{\text{PGA}} \\
\log_{10} D_n &= 0.22 - 10.12 \cdot r + 16.38 \cdot r^2 - 11.48 \cdot r^3 \\
\sigma_{\log_{10} D_n} &= 0.45
\end{aligned}
$$

##### Dn_JB07: Jibson (2007)

Dn_models.R:54-65:

$$
\begin{aligned}
\text{AI} &= \text{PGA}^{1.9228} \cdot \exp(2.6109) \quad \text{(m/s, if not provided)} \\
r &= \frac{k_y}{\text{PGA}} \\
\log_{10} D_n &= 0.561 \cdot \log_{10}(\text{AI}) - 3.833 \cdot \log_{10}(r) - 1.474 \\
\sigma_{\log_{10} D_n} &= 0.616
\end{aligned}
$$

##### Dn_SR08: Saygili & Rathje (2008)

Dn_models.R:71-83:

$$
\begin{aligned}
\text{AI} &= \text{PGA}^{1.9228} \cdot \exp(2.6109) \quad \text{(m/s, if not provided)} \\
r &= \frac{k_y}{\text{PGA}} \\
\mu_{\ln D_n} &= 2.39 - 5.24 \cdot r - 18.78 \cdot r^2 + 42.01 \cdot r^3 - 29.15 \cdot r^4 - 1.56 \cdot \ln(\text{PGA}) + 1.38 \cdot \ln(\text{AI}) \\
\sigma_{\ln D_n} &= 0.46 + 0.56 \cdot r
\end{aligned}
$$

Note: Heteroskedastic uncertainty (increases with $r$).

#### Flexible Sliding-Block Models

##### Dn_BT07: Bray & Travasarou (2007)

Dn_models.R:92-107:

$$
\begin{aligned}
\mu_{\ln D_n} &= -1.10 - 2.83 \cdot \ln(k_y) - 0.333 \cdot [\ln(k_y)]^2 + 0.566 \cdot \ln(k_y) \cdot \ln(\text{Sa}_{1.5T_s}) \\
&\quad + 3.04 \cdot \ln(\text{Sa}_{1.5T_s}) - 0.244 \cdot [\ln(\text{Sa}_{1.5T_s})]^2 + 0.278 \cdot (M_w - 7) \\
\sigma_{\ln D_n} &= 0.66
\end{aligned}
$$

where $\text{Sa}_{1.5T_s}$ is spectral acceleration at period $1.5 \cdot T_s$.

##### Dn_BM17: Bray & Macedo (2017)

Dn_models.R:116-142, for subduction/interface events:

$$
\begin{aligned}
a_0 &= c_1 + 0.550 \cdot M_w + c_2 \cdot T_s + c_3 \cdot T_s^2 - 3.353 \cdot \ln(k_y) - 0.390 \cdot [\ln(k_y)]^2 \\
a_1 &= 3.060 + 0.538 \cdot \ln(k_y) \\
a_2 &= -0.225 \\
\mu_{\ln D_n} &= a_0 + a_1 \cdot \ln(\text{Sa}_{1.5T_s}) + a_2 \cdot [\ln(\text{Sa}_{1.5T_s})]^2 \\
\sigma_{\ln D_n} &= 0.73
\end{aligned}
$$

Period-dependent coefficients:

$$
(c_1, c_2, c_3) = \begin{cases}
(-5.864, -9.421, 0.0) & \text{if } T_s < 0.1 \\
(-6.896, 3.081, -0.803) & \text{if } T_s \ge 0.1
\end{cases}
$$

##### Dn_BM19: Bray & Macedo (2019)

Dn_models.R:153-184, for shallow-crustal events with near-fault pulse effects:

$$
\begin{aligned}
I_{\text{pulse}} &= \begin{cases} 0 & \text{if PGV} \le 115 \text{ cm/s} \\ 1 & \text{if PGV} > 115 \text{ cm/s} \end{cases} \\
\Delta_{\text{PGV}} &= \begin{cases} 0 & \text{if PGV} \le 115 \\ -\ln(115) & \text{if PGV} > 115 \end{cases} \\
a_0 &= c_1 + 0.607 \cdot M_w + c_2 \cdot T_s + c_3 \cdot T_s^2 - 2.491 \cdot \ln(k_y) - 0.245 \cdot [\ln(k_y)]^2 \\
&\quad + I_{\text{pulse}} \cdot \ln(\text{PGV}) + \Delta_{\text{PGV}} \\
a_1 &= 2.703 + 0.344 \cdot \ln(k_y) \\
a_2 &= -0.089 \\
\mu_{\ln D_n} &= a_0 + a_1 \cdot \ln(\text{Sa}_{1.3T_s}) + a_2 \cdot [\ln(\text{Sa}_{1.3T_s})]^2 \\
\sigma_{\ln D_n} &= 0.74
\end{aligned}
$$

where $\text{Sa}_{1.3T_s}$ is spectral acceleration at period $1.3 \cdot T_s$.

Period-dependent coefficients:

$$
(c_1, c_2, c_3) = \begin{cases}
(-4.551, -9.690, 0.0) & \text{if } T_s < 0.1 \\
(-5.894, 3.152, -0.910) & \text{if } T_s \ge 0.1
\end{cases}
$$

PGV estimation when not provided (fitDn.R:58):

$$
\text{PGV} = \text{PGA}^{1.0529} \cdot \exp(0.1241) \cdot 100 \quad \text{[cm/s]}
$$

---

### Monte Carlo Sampling and Spectral Correlation

#### Required Spectral Ordinates

fitDn.R:31-35:

$$
\begin{aligned}
T_{d,0} &= 0.0 \quad \text{(PGA)} \\
T_{d,1.0} &= 1.0 \cdot T_s \\
T_{d,1.3} &= 1.3 \cdot T_s \\
T_{d,1.5} &= 1.5 \cdot T_s
\end{aligned}
$$

#### Baker & Jayaram (2009) Correlation

rhoBJ.R:25-36:

$$
\rho(T_n, T_0) = \exp\left(-\alpha \cdot \left|\ln\left(\frac{T_0}{T_n}\right)\right|\right)
$$

where $\alpha$ depends on rupture distance:

$$
\alpha = \begin{cases}
0.35 & \text{if } R_{\text{rup}} < 10 \text{ km} \\
0.30 & \text{if } 10 \le R_{\text{rup}} \le 20 \text{ km} \\
0.25 & \text{if } R_{\text{rup}} > 20 \text{ km}
\end{cases}
$$

#### Correlated Sampling Procedure

sampleSaCorr.R:15-50:

1. Build correlation matrix $\mathbf{C}$ with reference period (PGA) in position [1,1]:

$$
\mathbf{C}_{ij} = \begin{cases}
1 & \text{if } i = j \\
\rho(T_i, T_1) & \text{if } i = 1 \text{ or } j = 1 \\
0 & \text{otherwise (star-shaped)}
\end{cases}
$$

2. Ensure positive definiteness:

```r
ev <- eigen(C)$values
if (min(ev) <= 0) {
  C <- C + diag(abs(min(ev)) + 1e-8)
}
```

3. Sample Gaussian copula: $\mathbf{Z} \sim \mathcal{N}(\mathbf{0}, \mathbf{C})$

4. Transform to uniforms: $\mathbf{U} = \Phi(\mathbf{Z})$

5. Map to $\ln(\text{Sa})$ via empirical inverse quantile function $Q(u)$ (Hyman monotone spline, buildQSpline.R:47-51):

$$
\ln(\text{Sa}_i) = Q_i(U_i)
$$

6. Re-center to preserve tabulated mean:

$$
\text{Sa}_i \leftarrow \text{Sa}_i \cdot \frac{\mu_{\text{target}}}{\text{mean}(\text{Sa}_i)}
$$

#### Ensemble Displacement Sampling

fitDn.R:73-97:

For each model $m \in \{\text{YG91, AM88, JB07, BT07, SR08, BM17, BM19}\}$:

1. Evaluate $\mu_{\ln D_n,m}$ and $\sigma_{\ln D_n,m}$ from correlated $(PGA, Sa_{1.0T_s}, Sa_{1.3T_s}, Sa_{1.5T_s})$ samples
2. Draw $NS$ realizations: $D_{n,m,i} = \exp(\mathcal{N}(\mu_{\ln D_n,m}, \sigma_{\ln D_n,m}))$

Weighted quantiles (fitDn.R:120-127):

$$
D_n(p) = Q_{\text{weighted}}\left(\bigcup_{m} \{D_{n,m,i}\}, \; \text{weights} = \{w_m\}, \; \text{probs} = p\right)
$$

using Hmisc::wtd.quantile with type = "quantile".

---

## Dependencies

Core packages:
- `data.table`: Fast data manipulation and grouping operations
- `mvtnorm`: Multivariate normal sampling for Gaussian copula
- `Hmisc`: Weighted quantiles and means
- `stats`: Spline interpolation, random number generation
- `readxl`: User-supplied hazard data import
- `stringr`: String parsing for OpenQuake headers
- `future`, `future.apply`: Parallel processing

---

## References

### Site Response

Stewart, J.P. & Seyhan, E. (2017). Semi-empirical nonlinear site amplification model for global application. *Earthquake Spectra*, **33**(1), 87-110. https://doi.org/10.1193/100614eqs151m

### Rigid Sliding-Block Models

Ambraseys, N.N. & Menu, J.M. (1988). Earthquake-induced ground displacements. *Earthquake Engineering & Structural Dynamics*, **16**(7), 985-1006.

Yegian, M.K., Marciano, E.A., & Ghahraman, V.G. (1991). Earthquake-induced permanent deformations: probabilistic approach. *Journal of Geotechnical Engineering*, **117**(1), 35-50.

Jibson, R.W. (2007). Regression models for estimating coseismic landslide displacement. *Engineering Geology*, **91**(2-4), 209-218.

Saygili, G. & Rathje, E.M. (2008). Empirical predictive models for earthquake-induced sliding displacements of slopes. *Journal of Geotechnical and Geoenvironmental Engineering*, **134**(6), 790-803.

### Flexible Sliding-Block Models

Bray, J.D. & Travasarou, T. (2007). Simplified procedure for estimating earthquake-induced deviatoric slope displacements. *Journal of Geotechnical and Geoenvironmental Engineering*, **133**(4), 381-392.

Bray, J.D. & Macedo, J. (2017). Procedure for estimating shear-induced seismic slope displacement for shallow crustal earthquakes. *Journal of Geotechnical and Geoenvironmental Engineering*, **143**(12), 04017106.

Bray, J.D. & Macedo, J. (2019). 6th Ishihara lecture: Simplified procedure for estimating liquefaction-induced building settlement. *Soil Dynamics and Earthquake Engineering*, **102**, 215-231.

### Spectral Correlation

Baker, J.W. & Jayaram, N. (2009). Correlation of spectral acceleration values from NGA ground motion models. *Earthquake Spectra*, **24**(1), 299-317.

---

## Documentation

Function documentation is available via R help system:

```r
# Main functions
?buildGMDP    # Import hazard data
?fitSaF       # Site amplification
?fitDn        # Newmark displacements

# Individual displacement models
?Dn_AM88      # Ambraseys & Menu (1988)
?Dn_YG91      # Yegian et al. (1991)
?Dn_JB07      # Jibson (2007)
?Dn_SR08      # Saygili & Rathje (2008)
?Dn_BT07      # Bray & Travasarou (2007)
?Dn_BM17      # Bray & Macedo (2017)
?Dn_BM19      # Bray & Macedo (2019)

# Utilities
?sampleSaCorr        # Correlated spectral sampling
?rhoBJ               # Baker & Jayaram correlation
?buildQSpline        # Quantile function builder
?interpolateSaTable  # Period interpolation
```

For package overview:
```r
?newmark
```

---

## Contributing

Issues and pull requests are welcome at the [GitHub repository](https://github.com/averriK/newmark).

---

## License

See [LICENSE](LICENSE) file for details.

---

## Citation

When using this package in research, please cite:

```bibtex
@misc{newmark2024,
  author = {Verri Kozlowski, Alejandro},
  title = {newmark: Monte-Carlo Ensemble Newmark Displacements},
  year = {2024},
  version = {0.4.0},
  url = {https://github.com/averriK/newmark}
}
```

---

## Author

**Alejandro Verri Kozlowski**

- Email: averri@fi.uba.ar
- ORCID: [0000-0002-8535-1170](https://orcid.org/0000-0002-8535-1170)
- GitHub: [@averriK](https://github.com/averriK)

**Affiliation:**
- Facultad de Ingeniería, Universidad de Buenos Aires
