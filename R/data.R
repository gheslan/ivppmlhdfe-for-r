# data.R -- documentation for the bundled Monte-Carlo example datasets

#' Class A example data (individual + time fixed effects)
#'
#' A single Monte-Carlo realisation from the data-generating process of Kwon,
#' Larch, Yoon and Yotov (2026), Appendix C, for the Class A fixed-effect
#' structure \eqn{\psi_g = \alpha_i + \gamma_t}. The outcome is generated from
#' \eqn{y = \exp(\beta_1 x_1 + \beta_2 x_2 + \alpha_i + \gamma_t)\, v} with
#' \eqn{\beta_1 = 0.5}, \eqn{\beta_2 = 0.3}, and an endogenous regressor
#' \eqn{x_1} correlated with the multiplicative error \eqn{v}. \code{z} is a
#' valid excluded instrument for \code{x1}.
#'
#' @format A data frame with 1000 rows and 6 variables:
#' \describe{
#'   \item{id}{individual identifier (cross-sectional unit).}
#'   \item{year}{time identifier.}
#'   \item{z}{excluded instrument for \code{x1}.}
#'   \item{x2}{exogenous regressor.}
#'   \item{x1}{endogenous regressor (true coefficient 0.5).}
#'   \item{y}{non-negative count outcome.}
#' }
#' @source Kwon, Larch, Yoon and Yotov (2026), CESifo Working Paper 12641.
#' @examples
#' data(ivppmlhdfe_ClassA)
#' fit <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = ivppmlhdfe_ClassA)
#' coef(fit)["x1"]   # ~0.5, endogeneity removed
"ivppmlhdfe_ClassA"

#' Class B example data (two-way "gravity" fixed effects)
#'
#' A single Monte-Carlo realisation for the Class B fixed-effect structure
#' \eqn{\psi_g = \alpha_{it} + \gamma_{jt}} (exporter-time and importer-time
#' effects), as in Kwon, Larch, Yoon and Yotov (2026), Appendix C.
#'
#' @format A data frame with 8700 rows and 8 variables:
#' \describe{
#'   \item{pair}{directed exporter-importer pair identifier.}
#'   \item{year}{time identifier.}
#'   \item{exp}{exporter (origin) identifier.}
#'   \item{imp}{importer (destination) identifier.}
#'   \item{z}{excluded instrument for \code{x1}.}
#'   \item{x2}{exogenous regressor.}
#'   \item{x1}{endogenous regressor (true coefficient 0.5).}
#'   \item{y}{non-negative count outcome.}
#' }
#' @source Kwon, Larch, Yoon and Yotov (2026), CESifo Working Paper 12641.
#' @examples
#' data(ivppmlhdfe_ClassB)
#' fit <- ivppmlhdfe(y ~ x2 | exp^year + imp^year | x1 ~ z,
#'                   data = ivppmlhdfe_ClassB, cluster = ~pair)
"ivppmlhdfe_ClassB"

#' Class C example data (three-way "gravity" fixed effects)
#'
#' A single Monte-Carlo realisation for the Class C fixed-effect structure
#' \eqn{\psi_g = \alpha_{it} + \gamma_{jt} + \eta_{ij}} (exporter-time,
#' importer-time, and pair effects), as in Kwon, Larch, Yoon and Yotov (2026),
#' Appendix C.
#'
#' @format A data frame with 3800 rows and 8 variables:
#' \describe{
#'   \item{pair}{directed exporter-importer pair identifier.}
#'   \item{year}{time identifier.}
#'   \item{exp}{exporter (origin) identifier.}
#'   \item{imp}{importer (destination) identifier.}
#'   \item{z}{excluded instrument for \code{x1}.}
#'   \item{x2}{exogenous regressor.}
#'   \item{x1}{endogenous regressor (true coefficient 0.5).}
#'   \item{y}{non-negative count outcome.}
#' }
#' @source Kwon, Larch, Yoon and Yotov (2026), CESifo Working Paper 12641.
#' @examples
#' data(ivppmlhdfe_ClassC)
#' fit <- ivppmlhdfe(y ~ x2 | exp^year + imp^year + pair | x1 ~ z,
#'                   data = ivppmlhdfe_ClassC, cluster = ~pair)
"ivppmlhdfe_ClassC"
