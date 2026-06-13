# ivppmlhdfe-package.R -- package-level documentation and imports

#' ivppmlhdfe: IV-PPML with High-Dimensional Fixed Effects
#'
#' An R translation of the Stata command \code{ivppmlhdfe}. The package fits
#' instrumental-variable Poisson pseudo-maximum-likelihood (IV-PPML) models
#' with multiple sets of high-dimensional fixed effects, solving the additive
#' moment condition \eqn{E[q(y-\mu)] = 0} of Windmeijer and Santos Silva (1997)
#' via iteratively reweighted two-stage least squares (Correia, Guimaraes and
#' Zylkin 2020). It also provides the split-panel-jackknife bias correction of
#' Kwon, Larch, Yoon and Yotov (2026).
#'
#' @section Attribution:
#' This package is **only an R translation** of the Stata package `ivppmlhdfe`
#' by Ohyun Kwon, Mario Larch, Jangsu Yoon and Yoto V. Yotov, developed in their
#' working paper *Instrumental-Variable Poisson PML with High-Dimensional Fixed
#' Effects* (CESifo Working Paper 12641, 2026), available at
#' \url{https://www.ifo.de/en/cesifo/publications/2026/working-paper/instrumental-variable-poisson-pml-high-dimensional-fixed-effects}.
#'
#' All intellectual credit for the estimator, the methodology and the original
#' implementation belongs entirely to those authors. The R translation was
#' carried out by Guewen Heslan (PhD Candidate in Economics, Nantes Universite,
#' France) using Claude Opus, and is distributed at
#' \url{https://github.com/gheslan/ivppmlhdfe-for-r}; it adds no new methodology
#' and claims no credit for the underlying work.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{ivppmlhdfe}}}{formula interface (main entry point).}
#'   \item{\code{\link{ivppml_fit}}}{matrix interface (fast, for loops).}
#'   \item{\code{\link{ivppml_spj}}}{split-panel-jackknife bias correction.}
#'   \item{\code{\link{ivppml_tagsep}}}{separation tagging without estimation.}
#' }
#'
#' @section Numerical equivalence:
#' Plain PPML (no IV part) reproduces \code{ppmlhdfe} / \code{fixest::fepois}
#' point estimates to machine precision. The robust and cluster-robust sandwich
#' matrices are identical to those estimators up to the small-sample multiplier:
#' \code{ivppmlhdfe} uses \eqn{N/(N-1)} (robust) and \eqn{G/(G-1)} (cluster),
#' matching the Stata command exactly.
#'
#' @references
#' Kwon, O., M. Larch, J. Yoon, and Y.V. Yotov (2026). Instrumental-Variable
#' Poisson PML with High-Dimensional Fixed Effects. CESifo Working Paper 12641.
#'
#' Windmeijer, F.A.G. and J.M.C. Santos Silva (1997). Endogeneity in count data
#' models. \emph{Journal of Applied Econometrics} 12(3), 281--294.
#'
#' Correia, S., P. Guimaraes, and T. Zylkin (2020). Fast Poisson estimation with
#' high-dimensional fixed effects. \emph{Stata Journal} 20(1), 95--115.
#'
#' @keywords internal
#' @aliases ivppmlhdfe-package
#' @importFrom stats ave complete.cases model.frame model.matrix na.pass
#' @importFrom stats pchisq pnorm printCoefmat qnorm quantile sd setNames terms
#' @importFrom stats as.formula runif
#' @importFrom stats coef vcov nobs logLik fitted deviance residuals confint
#' @importFrom stats predict
#' @importFrom fixest demean
"_PACKAGE"
