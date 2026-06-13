# ivppmlhdfe

**Instrumental-Variable Poisson PML with High-Dimensional Fixed Effects** — an
R translation of the Stata command `ivppmlhdfe` by Kwon, Larch, Yoon and Yotov
(2026, CESifo Working Paper 12641).

> ### Attribution
>
> This package is **only a translation**. It is an R implementation of the Stata
> package by **Ohyun Kwon, Mario Larch, Jangsu Yoon, and Yoto V. Yotov**,
> developed in their working paper *Instrumental-Variable Poisson PML with
> High-Dimensional Fixed Effects*
> (<https://www.ifo.de/en/cesifo/publications/2026/working-paper/instrumental-variable-poisson-pml-high-dimensional-fixed-effects>).
>
> **All credit goes to the original authors.** The estimator, the methodology
> and the original implementation are entirely theirs; this package adds no new
> methodology.
>
> The R translation was carried out by **Guewen Heslan** (PhD Candidate in
> Economics, Nantes Université, France) using **Claude Opus**, and is available
> on GitHub at **[`gheslan/ivppmlhdfe-for-r`](https://github.com/gheslan/ivppmlhdfe-for-r)**.
>
> *Note: the importable R package is named `ivppmlhdfe` (R package names cannot
> contain hyphens); the GitHub repository is `ivppmlhdfe-for-r`.*


The estimator fits instrumental-variable Poisson pseudo-maximum-likelihood
(IV-PPML) models with multiple sets of high-dimensional fixed effects. It solves
the additive moment condition

```
E[q (y - mu)] = 0,   q = (x', z')',   mu = exp(x'b + d)
```

of Windmeijer and Santos Silva (1997) — replacing the endogenous regressor in
the PPML score with excluded instruments — via iteratively reweighted two-stage
least squares, with the fixed effects concentrated out at every iteration
following Correia, Guimaraes and Zylkin (2020). It also implements the
split-panel-jackknife (SPJ) bias correction with bootstrap standard errors that
the paper develops to remove the incidental-parameter bias that the IV moment
re-introduces.

## Installation

```r
# install.packages("remotes")
remotes::install_github("gheslan/ivppmlhdfe-for-r")
```

The only substantive dependency is **fixest**, used for the high-dimensional
fixed-effect demeaning (the same alternating-projections algorithm family as
Stata's `reghdfe`); everything else is base R.

## Quick start

```r
library(ivppmlhdfe)
data(ivppmlhdfe_ClassA)

# Class A: individual + time fixed effects, x1 endogenous, instrument z
m <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = ivppmlhdfe_ClassA)
summary(m)
coef(m)["x1"]      # ~0.5 (the true coefficient; endogeneity removed)
```

## Formula syntax

```
y ~ exog | fe1 + fe2 | endog ~ instruments
```

The two parts after the first `|` may appear in either order; the part
containing a `~` is the IV part, the other lists the fixed effects. Interacted
fixed effects use `^` (Stata's `#`).

| Specification | Formula |
|---|---|
| Class A (indiv + time) | `y ~ x2 \| id + year \| x1 ~ z` |
| Class B (two-way gravity) | `y ~ x2 \| exp^year + imp^year \| x1 ~ z` |
| Class C (three-way gravity) | `y ~ x2 \| exp^year + imp^year + pair \| x1 ~ z` |
| Plain PPML (no IV) | `y ~ x1 + x2 \| id + year` |
| IV, no fixed effects | `y ~ x2 \| x1 ~ z` |
| No exogenous regressors | `y ~ 1 \| id + year \| x1 ~ z` |

## Mapping of Stata options

| Stata | R |
|---|---|
| `absorb(...)` | FE part of the formula |
| `noabsorb` | omit the FE part |
| `vce(robust)` | default |
| `vce(cluster c1 c2)` | `cluster = ~c1 + c2` |
| `[pw=w]` / `[fw=w]` | `weights = ~w`, `weight_type =` |
| `exposure(v)` / `offset(v)` | `exposure = ~v` / `offset = ~v` |
| `tolerance(#)` / `itolerance(#)` | `tol` / `itol` |
| `maxiterations(#)` | `maxiter` |
| `separation(...)` | `separation =` |
| `guess()`, `standardize`, `keepsingletons` | same names |
| `tagsep(v)` | `ivppml_tagsep()` |
| `eform` / `irr` | `summary(fit, eform = TRUE)` |

## Cluster-robust and two-way clustering

```r
data(ivppmlhdfe_ClassB)
m <- ivppmlhdfe(y ~ x2 | exp^year + imp^year | x1 ~ z,
                data = ivppmlhdfe_ClassB, cluster = ~pair)        # one-way
m2 <- ivppmlhdfe(y ~ x2 | exp^year + imp^year | x1 ~ z,
                 data = ivppmlhdfe_ClassB, cluster = ~exp + imp)  # CGM two-way
```

## Split-panel jackknife bias correction

```r
s <- ivppml_spj(y ~ x2 | id + year | x1 ~ z, data = ivppmlhdfe_ClassA,
                class = "A", id = "id", time = "year",
                B = 1000, seed = 1)   # B = 0 for the point estimate only
s$b_spj        # bias-corrected estimate
s$se_ci_implied  # CI-implied bootstrap standard error (paper's aggregator)
```

The SPJ formulas follow the paper:

* **Class A** `3*b - mean(time halves) - mean(cross-section halves)`
* **Class B** `2*b - mean(4 country sub-panels)`
* **Class C** `4*b - 2*mean(country) - 2*mean(time) + mean(8 cells)`
* **Interacted single family** `2*b - mean(cross-section halves)`

## Numerical equivalence with Stata / ppmlhdfe

* **Point estimates.** Plain PPML (no IV part) reproduces `ppmlhdfe` /
  `fixest::fepois` coefficients to machine precision (verified to 1e-10–1e-7
  across the three fixed-effect classes). IV-PPML uses the identical IRLS engine
  with the instrument substituted into the slope moment.
* **Standard errors.** The robust and cluster-robust sandwich matrices (bread
  and meat) are identical to `ppmlhdfe`/`fepois`; the only difference is the
  small-sample multiplier, which `ivppmlhdfe` fixes at `N/(N-1)` (robust) and
  `G/(G-1)` (cluster), exactly matching the Stata command. With `fixest`'s
  small-sample correction switched off, the standard errors coincide to
  6+ digits.
* **Conventions reproduced.** singleton dropping, fixed-effect / ReLU
  separation, two-stage collinearity removal, absorbed-FE degrees of freedom
  (`e(df_a)` with the connected-components correction), the offset-centering
  trick, `mu` censoring, step-halving and adaptive inner tolerance all mirror
  `ivppmlhdfe.ado` / `ppmlhdfe`.

## References

Kwon, O., M. Larch, J. Yoon, and Y.V. Yotov (2026). *Instrumental-Variable
Poisson PML with High-Dimensional Fixed Effects.* CESifo Working Paper 12641.

Windmeijer, F.A.G. and J.M.C. Santos Silva (1997). Endogeneity in count data
models. *Journal of Applied Econometrics* 12(3), 281–294.

Correia, S., P. Guimaraes, and T. Zylkin (2020). Fast Poisson estimation with
high-dimensional fixed effects. *Stata Journal* 20(1), 95–115.

Weidner, M. and T. Zylkin (2021). Bias and consistency in three-way gravity
models. *Journal of International Economics* 132, 103513.
