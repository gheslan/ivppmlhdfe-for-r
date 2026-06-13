# methods.R -- S3 methods for "ivppmlhdfe" objects
#
# Display formatting follows the Replay program of ivppmlhdfe.ado; predict
# statistics follow ivppmlhdfe_p.ado.

#' @export
coef.ivppmlhdfe <- function(object, ...) object$coefficients

#' @export
vcov.ivppmlhdfe <- function(object, ...) object$vcov

#' @export
nobs.ivppmlhdfe <- function(object, ...) object$N

#' @export
logLik.ivppmlhdfe <- function(object, ...) {
  structure(object$ll, df = object$df_m + object$df_a,
            nobs = object$N, class = "logLik")
}

#' @export
fitted.ivppmlhdfe <- function(object, ...) object$fitted.values

#' @export
deviance.ivppmlhdfe <- function(object, ...) object$deviance

#' @export
residuals.ivppmlhdfe <- function(object,
                                 type = c("response", "pearson", "working",
                                          "deviance", "anscombe"),
                                 ...) {
  type <- match.arg(type)
  y <- object$y; mu <- object$fitted.values
  switch(type,
    response = y - mu,
    pearson  = (y - mu) / sqrt(mu),
    working  = (y - mu) / mu,
    anscombe = 1.5 * (y^(2/3) - mu^(2/3)) / mu^(1/6),
    deviance = {
      d <- 2 * ifelse(y > 0, mu - y + y * log(y / mu), mu)
      sign(y - mu) * sqrt(pmax(d, 0))
    })
}

#' @export
confint.ivppmlhdfe <- function(object, parm, level = 0.95, ...) {
  b <- coef(object)
  se <- sqrt(abs(diag(vcov(object))))
  if (missing(parm)) parm <- names(b)
  z <- stats::qnorm(1 - (1 - level) / 2)
  ci <- cbind(b - z * se, b + z * se)
  colnames(ci) <- paste(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2),
                               trim = TRUE), "%")
  rownames(ci) <- names(b)
  ci[parm, , drop = FALSE]
}

#' Predictions from an IV-PPML-HDFE fit
#'
#' Mirrors the `predict` options of the Stata command (`ivppmlhdfe_p.ado`).
#'
#' @param object an `"ivppmlhdfe"` fit.
#' @param type one of:
#'   \describe{
#'     \item{`"mu"` (default), `"n"`}{predicted mean `exp(xb + d [+ offset])`.}
#'     \item{`"xb"`}{linear prediction `X b + _cons` (+ offset if present).}
#'     \item{`"xbd"`, `"eta"`}{full linear predictor `xb + d`.}
#'     \item{`"d"`}{sum of the absorbed fixed effects.}
#'     \item{`"residuals"`, `"response"`, `"scores"`}{`y - mu`.}
#'     \item{`"pearson"`}{`(y - mu)/sqrt(mu)`.}
#'     \item{`"anscombe"`}{`1.5 (y^(2/3) - mu^(2/3)) / mu^(1/6)`.}
#'     \item{`"working"`}{`(y - mu)/mu`.}
#'     \item{`"deviance"`}{deviance contribution
#'       `2 (mu - y + y log(y/mu))` (and `2 mu` at `y = 0`).}
#'     \item{`"stdp"`}{standard error of the linear prediction from the
#'       slope-block VCE (the `_cons` variance is 0 by convention, so this
#'       reflects only the slope contribution).}
#'   }
#' @param ... unused.
#' @return numeric vector over the estimation sample.
#' @export
predict.ivppmlhdfe <- function(object,
                               type = c("mu", "n", "xb", "xbd", "eta", "d",
                                        "residuals", "response", "scores",
                                        "pearson", "anscombe", "working",
                                        "deviance", "stdp"),
                               ...) {
  type <- match.arg(type)
  if (type == "n") type <- "mu"
  if (type == "eta") type <- "xbd"

  y <- object$y
  mu <- object$fitted.values
  eta <- object$eta
  d <- object$d
  off <- if (is.null(object$offset)) 0 else object$offset
  b <- object$coefficients
  K <- length(b) - 1L

  xb <- eta - d                       # X b + _cons + offset
  switch(type,
    mu  = mu,
    xb  = xb,
    xbd = eta,
    d   = d,
    residuals = , response = , scores = y - mu,
    pearson  = (y - mu) / sqrt(mu),
    anscombe = 1.5 * (y^(2/3) - mu^(2/3)) / mu^(1/6),
    working  = (y - mu) / mu,
    deviance = 2 * ifelse(y > 0, mu - y + y * log(y / mu), mu),
    stdp = {
      X <- object$x
      if (is.null(X) || ncol(X) == 0L) {
        stop("type = \"stdp\" requires slope regressors in the model",
             call. = FALSE)
      }
      Vs <- object$vcov[seq_len(K), seq_len(K), drop = FALSE]
      sqrt(pmax(rowSums((X %*% Vs) * X), 0))
    })
}

#' @export
print.ivppmlhdfe <- function(x, eform = FALSE, digits = 4, ...) {
  print(summary(x, eform = eform), digits = digits, ...)
  invisible(x)
}

#' Summarize an IV-PPML-HDFE fit
#'
#' @param object an `"ivppmlhdfe"` fit.
#' @param eform,irr logical; display exponentiated coefficients
#'   (incidence-rate ratios). `irr` is an alias for `eform`.
#' @param ... unused.
#' @export
summary.ivppmlhdfe <- function(object, eform = FALSE, irr = FALSE, ...) {
  eform <- eform || irr
  b <- object$coefficients
  se <- sqrt(abs(diag(object$vcov)))
  zstat <- ifelse(se > 0, b / se, NA_real_)
  p <- 2 * stats::pnorm(-abs(zstat))
  lo <- b - stats::qnorm(0.975) * se
  hi <- b + stats::qnorm(0.975) * se
  if (eform) {
    tab <- cbind(`exp(b)` = exp(b), `Std. Err.` = exp(b) * se,
                 z = zstat, `P>|z|` = p,
                 `[95% lo` = exp(lo), `95% hi]` = exp(hi))
  } else {
    tab <- cbind(Coef. = b, `Std. Err.` = se, z = zstat, `P>|z|` = p,
                 `[95% lo` = lo, `95% hi]` = hi)
  }
  rownames(tab) <- names(b)
  out <- list(coeftable = tab, object = object, eform = eform)
  class(out) <- "summary.ivppmlhdfe"
  out
}

#' @export
print.summary.ivppmlhdfe <- function(x, digits = 4, ...) {
  o <- x$object
  K <- o$df_m
  cat(strrep("-", 78), "\n", sep = "")
  cat("IV-PPML with High-Dimensional Fixed Effects\n")
  cat(strrep("-", 78), "\n", sep = "")
  fmt <- function(lbl, val) cat(sprintf("%-22s %s\n", lbl, val))
  fmt("Dependent variable:", if (!is.null(o$depvar)) o$depvar else "y")
  if (length(o$endog_names)) {
    fmt("Endogenous:", paste(o$endog_names, collapse = " "))
    fmt("Instruments:", paste(o$inst_names, collapse = " "))
  }
  if (!is.null(o$absvars) && length(o$absvars)) {
    fmt("Absorbed FE:", paste(o$absvars, collapse = " "))
  }
  if (!is.null(o$offset_name)) fmt("Offset:", o$offset_name)
  cat(sprintf("No. of obs   = %10d    Residual df = %10d\n", o$N, o$df))
  cat(sprintf("Wald chi2(%d) = %10.2f    Prob > chi2 = %10.4f\n",
              K, o$chi2, stats::pchisq(o$chi2, K, lower.tail = FALSE)))
  cat(sprintf("Pseudo R2    = %10.4f    Deviance    = %10.4g\n",
              o$r2_p, o$deviance))
  cat(sprintf("Log pseudolikelihood = %.6g\n", o$ll))
  if (o$vcetype == "cluster") {
    cat(sprintf("No. of clusters = %d\n", o$N_clust))
  }
  cat(sprintf("Converged: %s  (iterations = %d)\n",
              if (o$converged) "yes" else "NO", o$iterations))
  if (o$num_singletons > 0) {
    cat(sprintf("(%d singleton observations dropped)\n", o$num_singletons))
  }
  if (o$num_separated > 0) {
    det <- c()
    if (o$num_sep_fe > 0) det <- c(det, paste0(o$num_sep_fe, " by fe"))
    if (o$num_sep_advanced > 0) det <- c(det, paste0(o$num_sep_advanced, " by simplex/relu"))
    if (o$num_sep_mu > 0) det <- c(det, paste0(o$num_sep_mu, " by mu"))
    cat(sprintf("(%d separated observations dropped: %s)\n",
                o$num_separated, paste(det, collapse = ", ")))
  }
  cat(strrep("-", 78), "\n", sep = "")
  stats::printCoefmat(x$coeftable, digits = digits, P.values = TRUE,
                      has.Pvalue = TRUE, cs.ind = 1:2, tst.ind = 3,
                      signif.stars = TRUE)
  cat(strrep("-", 78), "\n", sep = "")
  if (length(o$endog_names)) {
    cat("Endogenous: ", paste(o$endog_names, collapse = " "), "\n")
    cat("Instruments:", paste(o$inst_names, collapse = " "), "\n")
  }
  invisible(x)
}