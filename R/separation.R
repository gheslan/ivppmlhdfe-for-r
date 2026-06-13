# separation.R -- separation detection and singleton dropping
#
# Mirrors the Stata implementation:
#   * "fe" separation     : reghdfe iweight trick -- FE groups whose outcome
#                           is identically zero are dropped (interleaved with
#                           iterative singleton dropping).
#   * "relu" separation   : the ReLU / linear-programming-free certificate
#                           method of Correia, Guimaraes & Zylkin (2019),
#                           which detects separation caused by regressors
#                           and fixed effects jointly. This is the most
#                           general method; the "simplex" technique of
#                           ppmlhdfe is a fast special case whose detected
#                           set is contained in the ReLU set, so specifying
#                           "simplex" here is accepted and handled by the
#                           ReLU pass (the final estimation sample is
#                           identical).
#   * "mu" separation     : in-IRLS-loop detection of y == 0 observations
#                           whose linear predictor diverges to -Inf
#                           (implemented inside ivppml_fit()).

# ---------------------------------------------------------------------------
# Iterative singleton + zero-outcome-group dropping.
#
# Returns list(keep = logical(N), num_singletons, num_sep_fe).
# `drop_zero_groups = TRUE` reproduces reghdfe's iweight trick (weights = y):
# groups where sum(y) == 0 have zero total weight and are dropped.
# ---------------------------------------------------------------------------
fe_singleton_separation <- function(y, fe_ids, drop_singletons = TRUE,
                                    drop_zero_groups = TRUE) {
  N <- length(y)
  keep <- rep(TRUE, N)
  num_singletons <- 0L
  num_sep_fe <- 0L
  if (length(fe_ids) == 0L) {
    return(list(keep = keep, num_singletons = 0L, num_sep_fe = 0L))
  }
  repeat {
    changed <- FALSE
    for (fid in fe_ids) {
      idx <- which(keep)
      if (!length(idx)) break
      g <- group_id(fid[idx])
      if (drop_zero_groups) {
        gsum <- rowsum(y[idx], g)
        zero_g <- which(gsum[, 1L] == 0)
        if (length(zero_g)) {
          drop <- idx[g %in% zero_g]
          keep[drop] <- FALSE
          num_sep_fe <- num_sep_fe + length(drop)
          changed <- TRUE
          idx <- which(keep)
          if (!length(idx)) break
          g <- group_id(fid[idx])
        }
      }
      if (drop_singletons) {
        cnt <- tabulate(g)
        single <- which(cnt[g] == 1L)
        if (length(single)) {
          keep[idx[single]] <- FALSE
          num_singletons <- num_singletons + length(single)
          changed <- TRUE
        }
      }
    }
    if (!changed) break
  }
  list(keep = keep, num_singletons = num_singletons, num_sep_fe = num_sep_fe)
}

# ---------------------------------------------------------------------------
# ReLU separation (Correia, Guimaraes & Zylkin 2019, "Verifying the existence
# of maximum likelihood estimates for generalized linear models").
#
# Finds a certificate z* >= 0 with z*_i = 0 for y_i > 0 such that
# z* = X gamma + FE exactly. Observations with z*_i > 0 are separated.
#
# Implemented via the method-of-weighting: iterate weighted least squares of
# u on [X, FE] where y > 0 observations carry weight K >> 1 (forcing fitted
# values ~ 0 there), and y == 0 observations are rectified u <- max(uhat, 0).
#
# Returns list(separated = logical(N), certificate = numeric(N)).
# ---------------------------------------------------------------------------
relu_separation <- function(y, X, fe_df, w_user, tol = 1e-5, maxiter = 100L) {
  N <- length(y)
  is_zero <- y == 0
  if (!any(is_zero)) {
    return(list(separated = rep(FALSE, N), certificate = numeric(N)))
  }
  K <- if (is.null(X)) 0L else ncol(X)
  # method-of-weighting weights: large on y > 0
  Wbig <- N * 1e4
  ww <- ifelse(is_zero, 1, Wbig) * w_user
  u <- as.numeric(is_zero)
  uhat <- u
  for (it in seq_len(maxiter)) {
    M <- cbind(u, X)
    dm <- hdfe_demean(M, fe_df, weights = ww, tol = 1e-10)
    u_dm <- dm[, 1L]
    if (K > 0L) {
      X_dm <- dm[, -1L, drop = FALSE]
      XtWX <- crossprod(X_dm, ww * X_dm)
      XtWu <- crossprod(X_dm, ww * u_dm)
      g <- tryCatch(solve(XtWX, XtWu), error = function(e) {
        qr.coef(qr(sqrt(ww) * X_dm), sqrt(ww) * u_dm)
      })
      g[is.na(g)] <- 0
      r <- u_dm - drop(X_dm %*% g)
    } else {
      r <- u_dm
    }
    uhat <- u - r   # fitted value of the full (X + FE) regression
    # convergence: certificate found when uhat ~ u (i.e., u is exactly
    # representable) and uhat >= -tol everywhere
    if (max(abs(uhat - u)) < tol && min(uhat[is_zero]) > -tol) break
    # no separation if all fitted values on the zero sample are ~ 0
    if (max(uhat[is_zero]) < tol && max(abs(uhat[!is_zero])) < tol) {
      return(list(separated = rep(FALSE, N), certificate = numeric(N)))
    }
    u_new <- u
    u_new[is_zero] <- pmax(uhat[is_zero], 0)
    u_new[!is_zero] <- 0
    if (max(abs(u_new - u)) < tol * max(1, max(abs(u)))) { u <- u_new; break }
    u <- u_new
  }
  cert <- pmax(uhat, 0)
  cert[!is_zero] <- 0
  # rescale so that the largest certificate value is 1 (presentation only)
  mx <- max(cert)
  separated <- is_zero & cert > tol * max(1, mx)
  if (mx > 0) cert <- cert / mx
  list(separated = separated, certificate = cert)
}

# ---------------------------------------------------------------------------
# tagsep: run the full separation pipeline without estimating.
# ---------------------------------------------------------------------------

#' Tag separated observations without estimating
#'
#' Runs the same separation-detection pipeline as [ivppmlhdfe()] (singleton
#' dropping, fixed-effect separation, and ReLU separation) and returns the
#' tags without performing any estimation. Mirrors the `tagsep()` option of
#' the Stata command.
#'
#' @inheritParams ivppmlhdfe
#' @return A list with components:
#'   \item{separated}{logical vector (length `nrow(data)`); `TRUE` for
#'     observations that would be dropped as separated.}
#'   \item{singleton}{logical vector; `TRUE` for observations dropped as
#'     singletons.}
#'   \item{certificate}{numeric vector with the ReLU separation certificate
#'     (the `zvar` of the Stata command); positive entries identify
#'     separated observations.}
#'   \item{num_singletons, num_separated}{counts.}
#' @export
#' @examples
#' data(ivppmlhdfe_ClassA)
#' d <- ivppmlhdfe_ClassA
#' d$y[d$id == 1] <- 0  # make individual 1 all-zero (FE-separated)
#' tags <- ivppml_tagsep(y ~ x2 | id + year | x1 ~ z, data = d)
#' table(tags$separated)
ivppml_tagsep <- function(fml, data, weights = NULL,
                          separation = c("fe", "simplex", "relu"),
                          keep_singletons = FALSE) {
  parsed <- parse_ivppml_formula(fml, data)
  N <- nrow(data)
  keep0 <- parsed$complete
  y <- parsed$y[keep0]
  X <- parsed$X[keep0, , drop = FALSE]
  fe_ids <- lapply(parsed$fe_df, function(v) v[keep0])
  w_user <- if (is.null(weights)) rep(1, sum(keep0)) else {
    wv <- eval_fe_term(deparse(substitute(weights)), data)
    wv[keep0]
  }

  separation <- normalize_separation(separation)
  do_fe <- "fe" %in% separation
  do_relu <- any(c("relu", "simplex") %in% separation)

  res <- fe_singleton_separation(y, fe_ids,
                                 drop_singletons = !keep_singletons,
                                 drop_zero_groups = do_fe)
  keep <- res$keep
  cert <- numeric(length(y))
  sep_relu <- rep(FALSE, length(y))
  if (do_relu && any(keep)) {
    idx <- which(keep)
    fe_df_sub <- as.data.frame(lapply(fe_ids, function(v) v[idx]),
                               optional = TRUE)
    rl <- relu_separation(y[idx], X[idx, , drop = FALSE], fe_df_sub,
                          w_user[idx])
    sep_relu[idx] <- rl$separated
    cert[idx] <- rl$certificate
  }

  separated_full <- rep(FALSE, N)
  singleton_full <- rep(FALSE, N)
  cert_full <- rep(NA_real_, N)
  i0 <- which(keep0)
  # decompose drops: classify singleton vs separated by re-running
  res_sing_only <- fe_singleton_separation(y, fe_ids,
                                           drop_singletons = !keep_singletons,
                                           drop_zero_groups = FALSE)
  singleton_full[i0] <- !res_sing_only$keep
  separated_full[i0] <- (!keep & res_sing_only$keep) | sep_relu
  cert_full[i0] <- cert
  list(separated = separated_full,
       singleton = singleton_full,
       certificate = cert_full,
       num_singletons = sum(singleton_full),
       num_separated = sum(separated_full))
}

# Normalize the separation option (aliases match ppmlhdfe / ivppmlhdfe.ado)
normalize_separation <- function(separation) {
  if (is.null(separation)) return(character(0))
  separation <- tolower(separation)
  separation[separation == "ir"] <- "relu"
  if (length(separation) == 1L) {
    if (separation %in% c("", "def", "default", "standard", "on", "auto"))
      return(c("fe", "simplex", "relu"))
    if (separation %in% c("all", "full"))
      return(c("fe", "simplex", "relu", "mu"))
    if (separation %in% c("no", "off", "none"))
      return(character(0))
  }
  valid <- c("fe", "simplex", "relu", "mu")
  bad <- setdiff(separation, valid)
  if (length(bad)) {
    stop("separation: unknown technique(s) ", paste(bad, collapse = ", "),
         "\n  valid: fe simplex relu mu (or aliases: default standard all none)",
         call. = FALSE)
  }
  separation
}