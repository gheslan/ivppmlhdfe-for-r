# ivppmlhdfe 0.9.4

*Only an R translation* of the Stata package by Ohyun Kwon, Mario Larch, Jangsu
Yoon and Yoto V. Yotov (CESifo Working Paper 12641, 2026). All credit goes to
the original authors. Translation by Guewen Heslan (using Claude Opus),
available at https://github.com/gheslan/ivppmlhdfe-for-r.

* Initial R release: a faithful translation of the Stata command `ivppmlhdfe`.
* IRLS-IV engine with high-dimensional fixed-effect absorption (via `fixest`),
  heteroskedasticity-robust and (multi-way) cluster-robust standard errors,
  singleton / fixed-effect / ReLU separation, two-stage collinearity removal,
  and absorbed-FE degrees of freedom matching `reghdfe`/`ppmlhdfe`.
* Split-panel jackknife (`ivppml_spj`) bias correction with bootstrap standard
  errors for Classes A, B, C and the interacted single-family case.
* Plain PPML reproduces `ppmlhdfe` / `fixest::fepois` to machine precision.
