# Macroeconomic Analysis in R: A Time Series Modeling Project

A time series analysis of monthly Indian data on the exchange rate, inflation
and industrial production growth, using a Vector Autoregression (VAR).

This is a corrected version of an earlier draft. Two things were wrong in the
first version:

- The IIP growth cycle was shifted one month relative to the differenced
  exchange rate and inflation series, so the VAR was estimated on misaligned
  data.
- The `vars::causality()` results were presented as pairwise tests (for
  example "inflation -> exchange rate"), but the function performs a block
  test: it checks whether one variable helps predict all the other equations
  at once.

The numbers below come from the corrected alignment, and the tests are
described as block tests.

## 1. Data

Monthly observations for India, January 2014 to March 2020 (75 months). The
sample stops in March 2020 because the data file jumps from March to June 2020
and is not continuous afterwards.

| Variable | Description |
| --- | --- |
| monthly_exc_rate | Monthly average INR per USD |
| Inflation | CPI inflation, percent |
| Actual_IIP | IIP growth, year-on-year, percent |

The "output gap" in the earlier version was actually the HP-filtered cycle of
year-on-year IIP growth (lambda = 14400). It is called `iip_cycle` here. It is
a growth cycle, not the output gap.

## 2. Stationarity

ADF (up to 6 lags chosen by AIC) and KPSS, both with an intercept and no
trend:

| Series | ADF stat | 5% crit | KPSS stat | 5% crit | Read |
| --- | --- | --- | --- | --- | --- |
| Exchange rate | -0.99 | -2.88 | 1.45 | 0.46 | I(1) |
| Inflation | -3.86 | -2.88 | 0.65 | 0.46 | mixed |
| IIP cycle | -4.50 | -2.88 | 0.08 | 0.46 | I(0) |

Two comments:

- The exchange rate looks I(1): ADF cannot reject a unit root and KPSS rejects
  stationarity.
- For inflation the two tests disagree. ADF rejects the unit root while KPSS
  rejects stationarity. I keep the first difference, which is the conservative
  choice, but the disagreement is worth knowing about.

## 3. Cointegration

ARDL bounds test on `monthly_exc_rate ~ Inflation + iip_cycle`, case 3, with
exact p-values:

- Bounds F: 0.90, p = 0.877
- Bounds t: -0.79, p = 0.895

No evidence of a stable long-run relationship, so a VAR on the differenced
series (with the cycle in levels) is a reasonable next step.

## 4. The VAR

AIC, HQ, SC and FPE all select one lag. The model is a VAR(1) with a constant.
The roots of the characteristic polynomial are 0.283, 0.283 and 0.277, all
inside the unit circle, so the estimated system is stable.

Estimated coefficients (standard errors in brackets):

| | d_exc_rate | d_inflation | iip_cycle |
| --- | --- | --- | --- |
| d_exc_rate(-1) | 0.189 (0.120) | 0.096 (0.089) | 0.003 (0.453) |
| d_inflation(-1) | -0.279 (0.153) | 0.288 (0.114)* | 1.185 (0.577)* |
| iip_cycle(-1) | 0.070 (0.039) | -0.031 (0.029) | 0.090 (0.148) |

`*` marks p < 0.05. The two visible relationships are inflation's own
persistence and the inflation lag in the cycle equation.

## 5. Diagnostics

| Test | Result |
| --- | --- |
| Portmanteau (serial correlation) | p = 0.70, none detected |
| ARCH (multivariate) | p = 0.067, borderline |
| Jarque-Bera (normality) | p < 2.2e-16, rejected |

The normality rejection is driven by March 2020, the first pandemic month.
Dropping that single observation, the same test gives p = 0.31. That is why
the robust and split-sample checks below matter.

## 6. Block Granger tests

| Variable tested | F | p (OLS) | p (HC) |
| --- | --- | --- | --- |
| d_inflation | 3.25 | 0.041 | 0.068 |
| iip_cycle | 2.19 | 0.115 | 0.203 |
| d_exc_rate | 0.58 | 0.562 | 0.395 |

Each row asks whether the variable in question helps predict the other two
equations together. These are not pairwise tests.

Reading the table:

- Inflation has some predictive content for the block under OLS covariance
  (p = 0.041). It weakens with heteroskedasticity-consistent errors (p = 0.068)
  and disappears when March 2020 is dropped (p = 0.19). I would call this
  suggestive, not established.
- The cycle and the exchange rate do not help predict the rest of the system
  in any specification.

## 7. Conclusion

The system shows weak short-run dynamics. Inflation has mild predictive
content for the other variables, but the result depends on the covariance
estimator and on a single pandemic observation, so it should not be treated as
a firm finding. The exchange rate is close to a random walk in this sample,
and the IIP growth cycle adds little once it is included in the same system.
There is no long-run relationship among the levels.

## Files

- `R/VAR.R` — full script; runs from a clean clone
- `data/filename.csv` — monthly data
- `plots/` — figures
- `session_info.txt` — R and package versions from the last run

## How to run

```r
install.packages(c("readr", "lubridate", "ggplot2", "patchwork",
                   "vars", "urca", "ARDL", "mFilter",
                   "sandwich", "sessioninfo"))
```

Then:

```r
source("R/VAR.R")
```

## References

- Pesaran, M. H., Shin, Y., & Smith, R. J. (2001). Bounds testing approaches
  to the analysis of level relationships. *Journal of Applied Econometrics*,
  16(3), 289-326.
- Pfaff, B. (2008). VAR, SVAR and SVEC models: Implementation within R
  package vars. *Journal of Statistical Software*, 27(4), 1-32.
