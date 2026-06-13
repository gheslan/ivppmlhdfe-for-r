# ivppmlhdfe.R -- formula interface (main user entry point)

#' IV-PPML with high-dimensional fixed effects
#'
#' Estimates instrumental-variable Poisson pseudo-maximum-likelihood
#' (IV-PPML) models with multiple sets of high-dimensional fixed effects,
#' replicating the Stata command \command{ivppmlhdfe} (v0.9.4). The
#' estimator targets the additive moment condition
#' \eqn{E[q (y - \mu)] = 0} with \eqn{q = (x', z')'} (exogenous regressors
#' stacked with excluded instruments), together with the per-group
#' fixed-effect score \eqn{\sum_{g \in r} (y_g - \mu_g) = 0}, following
#' Windmeijer and Santos Silva (1997). It is solved via iteratively
#' reweighted 2SLS (Correia, Guimaraes and Zylkin 2020) with fixed effects
#' concentrated out at each iteration.
#'
#' @section Formula syntax:
#' The formula follows the multi-part convention
#' \preformatted{  y ~ exog | fe1 + fe2 | endog ~ instruments}
#' where the parts after the first `|` may appear in either order: the part
#' containing a `~` is the IV part, the other part lists the fixed effects
#' to absorb. Interacted fixed effects use `^` (e.g. `exp^year` is Stata's
#' `exp#year`). Examples:
#' \preformatted{
#'   y ~ x2 | id + year | x1 ~ z              # Class A
#'   y ~ x2 | exp^year + imp^year | x1 ~ z    # Class B
#'   y ~ x2 + x3 | id + year                  # plain PPML (no IV)
#'   y ~ x2 | x1 ~ z                          # IV, no FE (Stata noabsorb)
#'   y ~ 1  | id + year | x1 ~ z              # no exogenous regressors
#' }
#'
#' @section Mapping of Stata options:
#' \tabular{ll}{
#'   Stata \tab R \cr
#'   `absorb(...)` \tab FE part of `fml` \cr
#'   `noabsorb` \tab omit the FE part \cr
#'   `vce(robust)` \tab default \cr
#'   `vce(cluster c1 c2)` \tab `cluster = ~c1 + c2` \cr
#'   `[pw=w]` / `[fw=w]` \tab `weights = ~w`, `weight_type =` \cr
#'   `exposure(v)` / `offset(v)` \tab `exposure = ~v` / `offset = ~v` \cr
#'   `tolerance(#)` / `itolerance(#)` \tab `tol` / `itol` \cr
#'   `maxiterations(#)` \tab `maxiter` \cr
#'   `separation(...)` \tab `separation =` \cr
#'   `guess(...)`, `standardize`, `keepsingletons` \tab same names \cr
#'   `d(name)` / `d2` \tab always returned as `$d` \cr
#'   `tagsep(v)` \tab [ivppml_tagsep()] \cr
#'   `eform` / `irr` \tab `summary(fit, eform = TRUE)` \cr
#'   `verbose(#)` / `nolog` \tab `verbose` \cr
#' }
#'
#' @param fml a multi-part formula; see Details.
#' @param data a `data.frame` containing all variables.
#' @param weights optional weights: a one-sided formula (`~w`), a column
#'   name, or a numeric vector.
#' @param weight_type `"pweight"` (default; Stata `[pw=]`) or `"fweight"`
#'   (Stata `[fw=]`: frequency weights, effective N is `sum(weights)`).
#' @param cluster cluster specification for cluster-robust standard errors:
#'   a one-sided formula (`~pair`, `~exp + imp` for two-way CGM clustering,
#'   `~exp^imp` for an interaction cluster), a character vector of column
#'   names, or a vector / data.frame of ids. Default `NULL` gives
#'   heteroskedasticity-robust (sandwich) standard errors.
#' @param exposure optional exposure variable (its log enters as an offset
#'   with coefficient 1); one-sided formula, column name, or vector. Must be
#'   strictly positive.
#' @param offset optional offset variable (coefficient constrained to 1).
#'   Cannot be combined with `exposure`.
#' @param subset optional logical/integer vector selecting rows of `data`
#'   (Stata `if`/`in`).
#' @inheritParams ivppml_fit
#'
#' @return An object of class `"ivppmlhdfe"`, a list including:
#'   \item{coefficients}{named vector of slope coefficients plus `_cons`
#'     (the constant is recovered after partialling out the fixed effects,
#'     following the ivreg2 `partial()` convention; its variance is 0).}
#'   \item{vcov}{variance-covariance matrix (cons row/col are 0).}
#'   \item{fitted.values, eta, d, residuals}{fitted mean \eqn{\mu}, linear
#'     predictor, fixed-effect sum (Stata `d()`), and response residuals
#'     \eqn{y - \mu} on the estimation sample.}
#'   \item{N, N_full, df, df_a, df_m, chi2, ll, ll_0, r2_p, deviance}{fit
#'     statistics matching the Stata `e()` results of the same names
#'     (`N_full` is `e(N_full)`, `df` is `e(df)`, `r2_p` is `e(r2_p)`).}
#'   \item{converged, iterations}{IRLS convergence status.}
#'   \item{num_singletons, num_sep_fe, num_sep_advanced, num_sep_mu}{counts
#'     of dropped observations.}
#'   \item{esample}{logical vector over rows of `data` flagging the
#'     estimation sample (Stata `e(sample)`).}
#'
#' @references
#' Kwon, O., M. Larch, J. Yoon, and Y.V. Yotov (2026).
#' \emph{Instrumental-Variable Poisson PML with High-Dimensional Fixed
#' Effects.} CESifo Working Paper 12641.
#'
#' Windmeijer, F.A.G. and J.M.C. Santos Silva (1997). Endogeneity in count
#' data models: An application to demand for health care.
#' \emph{Journal of Applied Econometrics} 12(3), 281--294.
#'
#' Correia, S., P. Guimaraes, and T. Zylkin (2020). Fast Poisson estimation
#' with high-dimensional fixed effects. \emph{Stata Journal} 20(1), 95--115.
#'
#' @seealso [ivppml_fit()] for the matrix interface, [ivppml_tagsep()],
#'   [ivppml_spj()] for split-panel-jackknife bias correction,
#'   [predict.ivppmlhdfe()].
#'
#' @examples
#' # Class A: individual + time fixed effects
#' data(ivppmlhdfe_ClassA)
#' m <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = ivppmlhdfe_ClassA)
#' summary(m)
#'
#' # Class B: two-way gravity, clustered by pair
#' data(ivppmlhdfe_ClassB)
#' m2 <- ivppmlhdfe(y ~ x2 | exp^year + imp^year | x1 ~ z,
#'                  data = ivppmlhdfe_ClassB, cluster = ~pair)
#'
#' # Plain PPML (no IV): identical to ppmlhdfe / fixest::fepois
#' m0 <- ivppmlhdfe(y ~ x1 + x2 | id + year, data = ivppmlhdfe_ClassA)
#' @export
ivppmlhdfe <- function(fml, data,
                       weights = NULL,
                       weight_type = c("pweight", "fweight"),
                       cluster = NULL,
                       exposure = NULL, offset = NULL,
                       subset = NULL,
                       tol = 1e-8, itol = -1, maxiter = 10000L,
                       guess = c("default", "simple", "mean"),
                       standardize = FALSE,
                       separation = c("fe", "simplex", "relu"),
                       keep_singletons = FALSE,
                       verbose = -1L) {
  cl <- match.call()
  weight_type <- match.arg(weight_type)
  guess <- match.arg(guess)
  if (!is.data.frame(data)) stop("`data` must be a data.frame")
  data <- as.data.frame(data)

  if (!is.null(subset)) data <- data[subset, , drop = FALSE]
  N0 <- nrow(data)

  parsed <- parse_ivppml_formula(fml, data)

  # Resolve auxiliary variable specs
  get_aux <- function(spec, what) {
    if (is.null(spec)) return(NULL)
    v <- if (inherits(spec, "formula")) {
      tl <- attr(stats::terms(spec), "term.labels")
      if (length(tl) != 1L) stop(what, " must reference a single variable")
      eval_fe_term(tl, data)
    } else if (is.character(spec) && length(spec) == 1L) {
      eval_fe_term(spec, data)
    } else {
      as.numeric(spec)
    }
    if (length(v) != N0) stop(what, " has wrong length")
    v
  }

  w_vec <- get_aux(weights, "weights")
  expo <- get_aux(exposure, "exposure")
  offs <- get_aux(offset, "offset")
  if (!is.null(expo) && !is.null(offs)) {
    stop("cannot specify both exposure and offset")
  }
  if (!is.null(expo)) {
    if (any(expo <= 0, na.rm = TRUE)) stop("exposure must be greater than zero")
    offs <- log(expo)
  }

  clust_df <- parse_cluster(cluster, data, N0)

  # Sample marking (markout): complete cases on all model variables
  cc <- parsed$complete
  if (!is.null(w_vec)) cc <- cc & !is.na(w_vec)
  if (!is.null(offs)) cc <- cc & !is.na(offs)
  if (!is.null(clust_df)) cc <- cc & stats::complete.cases(clust_df)

  if (!any(cc)) stop("no observations")

  idx <- which(cc)
  fit <- ivppml_fit(
    y = parsed$y[idx],
    X_exog = parsed$X_exog[idx, , drop = FALSE],
    X_endog = parsed$X_endog[idx, , drop = FALSE],
    Z_excl = parsed$Z[idx, , drop = FALSE],
    fe = lapply(parsed$fe_df, function(v) v[idx]),
    weights = if (!is.null(w_vec)) w_vec[idx] else NULL,
    weight_type = weight_type,
    offset = if (!is.null(offs)) offs[idx] else NULL,
    cluster = if (!is.null(clust_df)) clust_df[idx, , drop = FALSE] else NULL,
    tol = tol, itol = itol, maxiter = maxiter,
    guess = guess, standardize = standardize,
    separation = separation, keep_singletons = keep_singletons,
    verbose = verbose,
    colnames_x = c(colnames(parsed$X_exog), colnames(parsed$X_endog)))

  # Map esample back to rows of `data`
  es <- rep(FALSE, N0)
  es[idx] <- fit$esample
  fit$esample <- es
  fit$call <- cl
  fit$fml <- fml
  fit$depvar <- parsed$depvar
  fit$absvars <- parsed$fe_names
  fit$offset_name <- if (!is.null(expo)) {
    paste0("ln(", aux_name(exposure), ")")
  } else if (!is.null(offs)) aux_name(offset) else NULL
  fit
}

aux_name <- function(spec) {
  if (is.null(spec)) return(NULL)
  if (inherits(spec, "formula")) {
    return(attr(stats::terms(spec), "term.labels")[1L])
  }
  if (is.character(spec)) return(spec[1L])
  deparse(substitute(spec))
}

# ---------------------------------------------------------------------------
# Formula parser
#
# Parts separated by `|`. The first part is `y ~ exog`. Among the remaining
# parts (0, 1 or 2): the part containing a `~` is the IV part
# `endog ~ instruments`; the other part lists fixed effects (with `^`
# interactions). Mirrors the detection logic of ParseIV in ivppmlhdfe.ado
# (`(` before the first `,` => IV mode).
# ---------------------------------------------------------------------------
parse_ivppml_formula <- function(fml, data) {
  if (!inherits(fml, "formula")) stop("`fml` must be a formula")
  txt <- paste(deparse(fml, width.cutoff = 500L), collapse = " ")
  parts <- strsplit(txt, "|", fixed = TRUE)[[1L]]
  parts <- trimws(parts)
  if (length(parts) > 3L) stop("formula has too many `|` parts")

  main <- stats::as.formula(parts[1L], env = environment(fml))
  depvar <- all.vars(main[[2L]])
  if (length(depvar) != 1L) stop("a single dependent variable is required")
  y <- eval(main[[2L]], envir = data, enclos = environment(fml))

  fe_part <- NULL; iv_part <- NULL
  for (p in parts[-1L]) {
    if (grepl("~", p, fixed = TRUE)) {
      if (!is.null(iv_part)) stop("only one IV part (`endog ~ inst`) allowed")
      iv_part <- p
    } else {
      if (!is.null(fe_part)) stop("only one fixed-effects part allowed")
      fe_part <- p
    }
  }

  # exogenous regressors via model.matrix (handles factors, interactions)
  rhs1 <- parts[1L]
  rhs1 <- sub("^[^~]*~", "", rhs1)
  rhs1 <- trimws(rhs1)
  N <- nrow(data)
  if (rhs1 %in% c("1", "0", "-1", "")) {
    X_exog <- matrix(0, N, 0L)
    cc_x <- rep(TRUE, N)
  } else {
    fx <- stats::as.formula(paste("~", rhs1), env = environment(fml))
    mf <- stats::model.frame(fx, data = data, na.action = stats::na.pass)
    X_exog <- stats::model.matrix(fx, mf)
    keep_cols <- colnames(X_exog) != "(Intercept)"
    X_exog <- X_exog[, keep_cols, drop = FALSE]
    cc_x <- stats::complete.cases(mf)
  }

  # IV part
  if (!is.null(iv_part)) {
    ivf <- stats::as.formula(iv_part, env = environment(fml))
    lhs_terms <- ivf[[2L]]
    fend <- stats::as.formula(paste("~", deparse(lhs_terms, width.cutoff = 500L)),
                              env = environment(fml))
    mfd <- stats::model.frame(fend, data = data, na.action = stats::na.pass)
    X_endog <- stats::model.matrix(fend, mfd)
    X_endog <- X_endog[, colnames(X_endog) != "(Intercept)", drop = FALSE]
    cc_d <- stats::complete.cases(mfd)

    finst <- stats::as.formula(paste("~", deparse(ivf[[3L]], width.cutoff = 500L)),
                               env = environment(fml))
    mfz <- stats::model.frame(finst, data = data, na.action = stats::na.pass)
    Z <- stats::model.matrix(finst, mfz)
    Z <- Z[, colnames(Z) != "(Intercept)", drop = FALSE]
    cc_z <- stats::complete.cases(mfz)
  } else {
    X_endog <- matrix(0, N, 0L)
    Z <- matrix(0, N, 0L)
    cc_d <- cc_z <- rep(TRUE, N)
  }

  # FE part
  if (!is.null(fe_part) && !fe_part %in% c("0", "1", "")) {
    fe_terms <- strsplit(fe_part, "+", fixed = TRUE)[[1L]]
    fe_terms <- trimws(fe_terms)
    fe_terms <- fe_terms[nzchar(fe_terms)]
    fe_df <- lapply(fe_terms, function(tt) eval_fe_term(tt, data))
    names(fe_df) <- fe_terms
    cc_fe <- Reduce(`&`, lapply(fe_df, function(v) !is.na(v)))
  } else {
    fe_df <- list()
    fe_terms <- character(0)
    cc_fe <- rep(TRUE, N)
  }

  complete <- !is.na(y) & cc_x & cc_d & cc_z & cc_fe

  list(y = as.numeric(y), X_exog = X_exog, X_endog = X_endog, Z = Z,
       fe_df = fe_df, fe_names = fe_terms, depvar = depvar,
       complete = complete)
}