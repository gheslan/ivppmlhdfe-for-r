# fit.R -- low-level IRLS-IV engine
#
# This is a line-by-line faithful translation of the IRLS-IV loop in
# ivppmlhdfe.ado (Mata function ivppmlhdfe_irls) / the Julia reference
# implementation (IVPPMLFixedEffectModels.jl, ivppml_reg).
#
# Algorithm: iteratively reweighted 2SLS targeting the additive moment
# E[q (y - mu)] = 0, q = (x', z')' (Windmeijer & Santos Silva 1997), with
# fixed effects concentrated out at each iteration via weighted
# within-transformation (Correia, Guimaraes & Zylkin 2020).

#' Low-level IV-PPML-HDFE fitting engine
#'
#' Matrix-interface estimator used internally by [ivppmlhdfe()]. Most users
#' should call [ivppmlhdfe()]; this function is exported for programmatic
#' use (Monte Carlo loops, bootstrap) where formula parsing overhead matters.
#'
#' @param y numeric vector of non-negative outcomes.
#' @param X_exog numeric matrix of exogenous regressors (may have 0 columns).
#' @param X_endog numeric matrix of endogenous regressors (may have 0
#'   columns; if so the estimator degenerates to standard PPML).
#' @param Z_excl numeric matrix of excluded instruments.
#' @param fe list (or data.frame) of fixed-effect id vectors; empty list or
#'   `NULL` for an intercept-only model (Stata's `noabsorb`).
#' @param weights optional vector of weights (see `weight_type`).
#' @param weight_type `"pweight"` (default, sampling weights) or
#'   `"fweight"` (frequency weights). With frequency weights the effective
#'   number of observations is `sum(weights)` and the sandwich meat scales
#'   with one power of the weight, matching `ppmlhdfe`.
#' @param offset optional offset vector (coefficient constrained to 1).
#' @param cluster optional data.frame / list / vector of cluster ids for
#'   (multi-way) cluster-robust standard errors.
#' @param tol outer IRLS convergence tolerance (Stata `tolerance()`,
#'   default `1e-8`).
#' @param itol inner fixed-effect-solver target tolerance (Stata
#'   `itolerance()`); `-1` (default) uses `max(1e-12, 0.1 * tol)`.
#' @param maxiter maximum IRLS iterations (default 10000).
#' @param guess initial value rule for `mu`: `"default"`/`"simple"` uses
#'   `0.5 * (y + mean(y))`, `"mean"` uses `mean(y)`.
#' @param standardize logical; divide columns of X and Z by their weighted
#'   standard deviation during estimation for numerical stability
#'   (coefficients and VCE are back-transformed; results are identical up
#'   to numerical precision).
#' @param separation character vector of separation-detection techniques;
#'   any of `"fe"`, `"simplex"`, `"relu"`, `"mu"`, or the aliases
#'   `"default"` (= fe simplex relu), `"all"` (adds mu), `"none"`. The
#'   `"simplex"` technique is accepted for syntax compatibility and is
#'   subsumed by the `"relu"` pass (its detected set is a subset of the
#'   ReLU set, so the estimation sample is identical).
#' @param keep_singletons logical; keep singleton observations (default
#'   drops them iteratively, matching reghdfe).
#' @param verbose integer; `-1` silent (default), `0` iteration log,
#'   `1`+ adds detail.
#' @param colnames_x optional column names for `[X_exog, X_endog]` used in
#'   the returned coefficient vector.
#'
#' @return An object of class `"ivppmlhdfe"`; see [ivppmlhdfe()].
#' @references
#' Windmeijer, F.A.G. and J.M.C. Santos Silva (1997). Endogeneity in count
#' data models. \emph{Journal of Applied Econometrics} 12(3), 281--294.
#'
#' Correia, S., P. Guimaraes, and T. Zylkin (2020). Fast Poisson estimation
#' with high-dimensional fixed effects. \emph{Stata Journal} 20(1), 95--115.
#'
#' Kwon, O., M. Larch, J. Yoon, and Y.V. Yotov (2026).
#' Instrumental-Variable Poisson PML with High-Dimensional Fixed Effects.
#' CESifo Working Paper 12641.
#' @export
ivppml_fit <- function(y, X_exog = NULL, X_endog = NULL, Z_excl = NULL,
                       fe = NULL, weights = NULL,
                       weight_type = c("pweight", "fweight"),
                       offset = NULL, cluster = NULL,
                       tol = 1e-8, itol = -1, maxiter = 10000L,
                       guess = c("default", "simple", "mean"),
                       standardize = FALSE,
                       separation = c("fe", "simplex", "relu"),
                       keep_singletons = FALSE,
                       verbose = -1L,
                       colnames_x = NULL) {

  guess <- match.arg(guess)
  weight_type <- match.arg(weight_type)
  separation <- normalize_separation(separation)
  do_sep_fe <- "fe" %in% separation
  do_sep_relu <- any(c("relu", "simplex") %in% separation)
  do_sep_mu <- "mu" %in% separation

  # ---- Assemble and validate inputs -------------------------------------
  y <- as.numeric(y)
  N_input <- length(y)
  as_mat <- function(M) {
    if (is.null(M)) return(matrix(0, N_input, 0L))
    M <- as.matrix(M)
    storage.mode(M) <- "double"
    M
  }
  X_exog <- as_mat(X_exog); X_endog <- as_mat(X_endog); Z_excl <- as_mat(Z_excl)
  K_exog <- ncol(X_exog); K_endo <- ncol(X_endog); n_inst <- ncol(Z_excl)
  X <- cbind(X_exog, X_endog)
  K <- K_exog + K_endo

  if (is.null(colnames_x)) {
    cn_e <- colnames(X_exog); cn_d <- colnames(X_endog)
    if (is.null(cn_e) && K_exog > 0) cn_e <- paste0("x", seq_len(K_exog))
    if (is.null(cn_d) && K_endo > 0) cn_d <- paste0("d", seq_len(K_endo))
    colnames_x <- c(cn_e, cn_d)
  }
  inst_names <- colnames(Z_excl)
  if (is.null(inst_names) && n_inst > 0) inst_names <- paste0("z", seq_len(n_inst))

  if (any(!is.finite(y))) stop("depvar contains NA/NaN/Inf")
  if (any(y < 0)) stop("depvar must be nonnegative")
  if (sum(y) == 0) stop("depvar is identically zero; model is degenerate")
  if (K == 0L) {
    stop("must specify at least one regressor (an empty model is not supported)")
  }
  if (n_inst < K_endo) stop("equation not identified")
  if (K > 0 && any(!is.finite(X))) stop("regressors contain NA/NaN/Inf")
  if (n_inst > 0 && any(!is.finite(Z_excl))) stop("instruments contain NA/NaN/Inf")

  # Reject instruments identical to a regressor (degenerate first stage)
  if (n_inst > 0 && K > 0) {
    for (j in seq_len(n_inst)) {
      for (k in seq_len(K)) {
        if (identical(Z_excl[, j], X[, k])) {
          stop("instrument cannot equal a regressor: ",
               inst_names[j], call. = FALSE)
        }
      }
    }
  }
  # Reject constant (zero-variance) instruments
  if (n_inst > 0) {
    for (j in seq_len(n_inst)) {
      if (min(Z_excl[, j]) == max(Z_excl[, j])) {
        stop("instrument ", inst_names[j],
             " has zero variance (constant column)", call. = FALSE)
      }
    }
  }

  # Weights
  has_weight <- !is.null(weights)
  w_user <- if (has_weight) as.numeric(weights) else rep(1, N_input)
  if (any(!is.finite(w_user))) stop("weights contain NA/NaN/Inf")
  if (any(w_user < 0)) stop("negative weights encountered")

  # Offset
  has_offset <- !is.null(offset)
  offset_orig <- if (has_offset) as.numeric(offset) else rep(0, N_input)
  if (has_offset && any(!is.finite(offset_orig))) stop("offset contains NA/NaN/Inf")

  # FE ids
  if (is.null(fe)) fe <- list()
  if (is.data.frame(fe)) fe <- as.list(fe)
  fe_ids <- lapply(fe, function(v) {
    if (length(v) != N_input) stop("FE vector length mismatch")
    group_id(v)
  })
  n_fe <- length(fe_ids)

  # Cluster ids
  clust_df <- parse_cluster(cluster, data = NULL, n = N_input)
  has_cluster <- !is.null(clust_df) && ncol(clust_df) > 0

  # ---- Sample restriction ------------------------------------------------
  # Drop zero-weight rows (markout semantics for [pw=] / [fw=])
  keep <- rep(TRUE, N_input)
  if (has_weight) keep <- keep & (w_user > 0)

  # FE separation + singleton dropping (the reghdfe iweight trick)
  num_singletons <- 0L; num_sep_fe <- 0L
  if (n_fe > 0L) {
    idx <- which(keep)
    fs <- fe_singleton_separation(
      y[idx], lapply(fe_ids, function(v) v[idx]),
      drop_singletons = !keep_singletons,
      drop_zero_groups = do_sep_fe)
    keep[idx[!fs$keep]] <- FALSE
    num_singletons <- fs$num_singletons
    num_sep_fe <- fs$num_sep_fe
  }
  if (keep_singletons && verbose > -1L) {
    warning("keeping singleton groups will keep fixed effects that cause separation",
            call. = FALSE)
  }

  subset_all <- function(keep_now) {
    y <<- y[keep_now]
    X <<- X[keep_now, , drop = FALSE]
    Z_excl <<- Z_excl[keep_now, , drop = FALSE]
    w_user <<- w_user[keep_now]
    offset_orig <<- offset_orig[keep_now]
    fe_ids <<- lapply(fe_ids, function(v) v[keep_now])
    if (has_cluster) clust_df <<- clust_df[keep_now, , drop = FALSE]
  }
  esample <- keep
  subset_all(keep)
  N <- length(y)
  if (N == 0L) stop("no observations")

  # ---- ReLU / simplex separation (depends on y, X, FE -- not Z) ----------
  num_sep_advanced <- 0L
  if (do_sep_relu && K > 0L && any(y == 0)) {
    fe_df_now <- as.data.frame(fe_ids, optional = TRUE)
    rl <- relu_separation(y, X, if (n_fe) fe_df_now else NULL, w_user)
    if (any(rl$separated)) {
      keep2 <- !rl$separated
      num_sep_advanced <- sum(rl$separated)
      idx_es <- which(esample)
      esample[idx_es[!keep2]] <- FALSE
      subset_all(keep2)
      N <- length(y)
      # cascade: dropping separated obs can create new singletons /
      # zero-outcome groups
      if (n_fe > 0L) {
        fs2 <- fe_singleton_separation(y, fe_ids,
                                       drop_singletons = !keep_singletons,
                                       drop_zero_groups = do_sep_fe)
        if (!all(fs2$keep)) {
          idx_es <- which(esample)
          esample[idx_es[!fs2$keep]] <- FALSE
          subset_all(fs2$keep)
          N <- length(y)
          num_singletons <- num_singletons + fs2$num_singletons
          num_sep_fe <- num_sep_fe + fs2$num_sep_fe
        }
      }
    }
  }
  if (N == 0L) stop("no observations remaining after separation detection")

  # Re-code FE ids on the final sample and build the demeaning frame
  fe_ids <- lapply(fe_ids, group_id)
  fe_df <- if (n_fe) as.data.frame(fe_ids, optional = TRUE) else NULL
  if (!is.null(fe_df)) names(fe_df) <- paste0("fe", seq_len(n_fe))

  # G = 1 demotion: single cluster variable with a single group
  if (has_cluster && ncol(clust_df) == 1L) {
    if (length(unique(clust_df[[1L]])) <= 1L) {
      warning("only 1 cluster; switching to robust VCE", call. = FALSE)
      has_cluster <- FALSE
      clust_df <- NULL
    }
  }

  # Center offset to weighted mean zero for the IRLS loop (numerical
  # stability with large absolute-scale offsets; the fitted intercept
  # absorbs the shift). Matches ivppmlhdfe.ado and ppmlhdfe.mata:461.
  offset_mean <- if (has_offset) sum(w_user * offset_orig) / sum(w_user) else 0
  off <- offset_orig - offset_mean

  # ---- Adaptive tolerances -----------------------------------------------
  start_inner_tol <- 1e-4
  target_inner_tol <- if (itol > 0) itol else max(1e-12, 0.1 * tol)
  inner_tol <- max(start_inner_tol, tol)
  alt_tol <- start_inner_tol

  # ---- Two-stage collinearity removal ------------------------------------
  # Per-block (exog / endog / inst) rank detection AFTER projecting through
  # the FE structure with USER weights (not the IRLS weights).
  n_dropped_exog <- 0L; n_dropped_endog <- 0L; n_dropped_inst <- 0L
  if (K + n_inst > 0L) {
    data_raw <- cbind(X, Z_excl)
    data_coll <- hdfe_demean(data_raw, fe_df, weights = w_user,
                             tol = target_inner_tol)
    exog_keep <- if (K_exog > 0)
      select_not_collinear(data_coll[, seq_len(K_exog), drop = FALSE],
                           data_raw[, seq_len(K_exog), drop = FALSE], w_user)
      else integer(0)
    endog_keep <- if (K_endo > 0)
      select_not_collinear(data_coll[, K_exog + seq_len(K_endo), drop = FALSE],
                           data_raw[, K_exog + seq_len(K_endo), drop = FALSE],
                           w_user)
      else integer(0)
    inst_keep <- if (n_inst > 0)
      select_not_collinear(data_coll[, K + seq_len(n_inst), drop = FALSE],
                           data_raw[, K + seq_len(n_inst), drop = FALSE],
                           w_user)
      else integer(0)

    n_dropped_exog <- K_exog - length(exog_keep)
    n_dropped_endog <- K_endo - length(endog_keep)
    n_dropped_inst <- n_inst - length(inst_keep)

    if (n_dropped_exog + n_dropped_endog + n_dropped_inst > 0L) {
      x_idx <- c(exog_keep, K_exog + endog_keep)
      X <- X[, x_idx, drop = FALSE]
      Z_excl <- Z_excl[, inst_keep, drop = FALSE]
      colnames_x <- colnames_x[x_idx]
      inst_names <- inst_names[inst_keep]
      K_exog <- length(exog_keep); K_endo <- length(endog_keep)
      K <- K_exog + K_endo; n_inst <- length(inst_keep)
      if (n_inst < K_endo) {
        stop("equation not identified after collinearity removal (n_inst=",
             n_inst, " < K_endo=", K_endo, ")", call. = FALSE)
      }
      if (K == 0L) stop("all regressors dropped by collinearity removal")
      if (verbose > -1L) {
        if (n_dropped_exog > 0)
          message("note: ", n_dropped_exog,
                  " exog variable(s) omitted because of collinearity")
        if (n_dropped_endog > 0)
          message("note: ", n_dropped_endog,
                  " endog variable(s) omitted because of collinearity")
        if (n_dropped_inst > 0)
          message("note: ", n_dropped_inst,
                  " instrument(s) omitted because of collinearity")
      }
    }
  }
  L <- K_exog + n_inst

  # ---- Standardization (optional) ----------------------------------------
  stdev_x <- rep(1, K); stdev_z <- rep(1, n_inst)
  if (standardize && K > 0L) {
    sx <- wcolsd(X, w_user)
    stdev_x <- ifelse(sx > 0, sx, 1)
    X <- sweep(X, 2L, stdev_x, "/")
    if (n_inst > 0L) {
      sz <- wcolsd(Z_excl, w_user)
      stdev_z <- ifelse(sz > 0, sz, 1)
      Z_excl <- sweep(Z_excl, 2L, stdev_z, "/")
    }
  }

  # ---- Initialise mu ------------------------------------------------------
  mean_y <- sum(w_user * y) / sum(w_user)
  mu <- if (guess == "mean") rep(mean_y, N) else 0.5 * (y + mean_y)
  # y-dependent censoring (ppmlhdfe censor_mu): floor at max(0.05*y, 1e-3)
  mu <- pmax(mu, pmax(0.05 * y, 1e-3))
  eta <- log(mu) + off

  # ---- IRLS-IV loop --------------------------------------------------------
  converged <- FALSE; ok <- 0L
  deviance <- Inf; eps_val <- Inf; beta_change <- Inf
  b <- rep(0, K); b_old <- b
  max_step_halving <- 2L; step_memory <- 0.9; num_sh <- 0L
  sep_mask <- rep(FALSE, N); num_sep_mu <- 0L
  eps100 <- 100 * .Machine$double.eps
  iter <- 0L
  resid <- numeric(N)
  zero_idx <- which(y == 0)
  Z_dm <- NULL; X_dm <- NULL; z_dm <- NULL

  if (verbose > -1L) {
    cat(sprintf("\nIRLS-IV iterations (N = %d, K = %d, L = %d)\n", N, K, L))
    cat(strrep("-", 60), "\n", sep = "")
  }

  for (it in seq_len(maxiter)) {
    iter <- it

    # (a0) Runaway-divergence guard (Stata rc=9003)
    if (any(!is.finite(mu)) || any(!is.finite(eta))) {
      stop("mu/eta has infinite or missing values at iteration ", it,
           "; IRLS divergence. Try separation = \"all\", standardize = TRUE, ",
           "or a simpler FE structure.", call. = FALSE)
    }
    if (it > 10L && K > 0L) {
      if (max(abs(b)) > 1e6 || max(abs(eta)) > 30) {
        stop(sprintf(paste0(
          "IRLS diverged at iter %d (max|b|=%.2e, max|eta|=%.2e); aborting. ",
          "Try separation = \"all\", standardize = TRUE, or simpler FE."),
          it, max(abs(b)), max(abs(eta))), call. = FALSE)
      }
    }

    # (a) Working dependent variable
    z <- eta - off - 1 + y / mu
    if (do_sep_mu && length(zero_idx)) {
      z[zero_idx] <- eta[zero_idx] - off[zero_idx] - 1
    }

    # (b) IRLS weights (floored)
    irls_w <- pmax(w_user * mu, 1e-20)

    # (d-e) Demean [z, X, Z_excl]
    M <- cbind(z, X, Z_excl)
    M <- edit_to_zero(M, min(tol, 1e-12))
    M <- hdfe_demean(M, fe_df, weights = irls_w, tol = inner_tol)
    M <- edit_to_zero(M, min(tol, 1e-12))

    z_dm <- M[, 1L]
    X_dm <- M[, 1L + seq_len(K), drop = FALSE]
    Z_dm <- if (K_exog > 0L && n_inst > 0L) {
      cbind(M[, 1L + seq_len(K_exog), drop = FALSE],
            M[, 1L + K + seq_len(n_inst), drop = FALSE])
    } else if (K_exog > 0L) {
      M[, 1L + seq_len(K_exog), drop = FALSE]
    } else {
      M[, 1L + K + seq_len(n_inst), drop = FALSE]
    }

    # (g) Weighted 2SLS
    b_old <- b
    ZwZ <- crossprod(Z_dm, irls_w * Z_dm)
    ZwX <- crossprod(Z_dm, irls_w * X_dm)
    Pi <- tryCatch(solve(ZwZ, ZwX), error = function(e)
      stop("2SLS solve failed at iter ", it,
           " -- the design appears rank-deficient. Check for perfectly ",
           "collinear regressors/instruments or set standardize = TRUE.",
           call. = FALSE))
    Xhat <- Z_dm %*% Pi
    XhwX <- crossprod(Xhat, irls_w * X_dm)
    Xhwy <- crossprod(Xhat, irls_w * z_dm)
    b <- drop(tryCatch(solve(XhwX, Xhwy), error = function(e)
      stop("2SLS solve failed at iter ", it,
           " -- the design appears rank-deficient.", call. = FALSE)))
    resid <- z_dm - drop(X_dm %*% b)
    beta_change <- max(abs(b - b_old))

    # (h) Update eta
    old_eta <- eta
    eta <- z - resid + off

    # (h2) Mu-separation check (accumulating; once separated, always)
    if (do_sep_mu && it > 1L) {
      log_septol <- log(1e-6)
      pos <- y > 0
      adjusted <- if (any(pos)) log_septol + min(min(eta[pos]) + 5, 0)
                  else log_septol
      new_sep <- !sep_mask & y == 0 & eta <= adjusted
      if (any(new_sep)) {
        sep_mask <- sep_mask | new_sep
        n_new <- sum(new_sep)
        num_sep_mu <- sum(sep_mask)
        if (verbose > -1L)
          cat(sprintf("(mu-separation: %d obs detected at iter %d)\n", n_new, it))
      }
    }

    # (i) Update mu
    mu <- exp(eta)
    if (do_sep_mu && num_sep_mu > 0L) mu[sep_mask] <- 0
    mu <- pmax(mu, eps100)

    # (j) Deviance with ppmlhdfe numerical safeguards
    old_deviance <- deviance
    dev1 <- sum((mu - y) * w_user)
    ypos <- y > 0
    dev2 <- sum(w_user[ypos] * y[ypos] * (log(y[ypos]) - eta[ypos]))
    deviance <- dev1 + dev2
    if (2 * deviance / N < .Machine$double.eps) deviance <- 0
    deviance <- 2 * deviance
    if (abs(deviance) < .Machine$double.eps) deviance <- 0
    if (deviance < 0) deviance <- 0

    # (k) Convergence + step-halving + adaptive tolerance
    is_sh <- FALSE
    if (it > 1L) {
      delta_dev <- old_deviance - deviance
      if (deviance < 0.1 * delta_dev) delta_dev <- deviance
      denom <- max(min(deviance, old_deviance), 0.1)
      eps_val <- abs(delta_dev) / denom
      if (eps_val < tol) {
        if (inner_tol <= 1.1 * target_inner_tol || n_fe <= 1L) {
          ok <- ok + 1L
          if (ok >= 1L) converged <- TRUE
        }
      } else if (delta_dev < 0 && num_sh < max_step_halving) {
        eta <- step_memory * old_eta + (1 - step_memory) * eta
        if (num_sh > 0L) eta <- pmax(eta, -10)   # clip (ppmlhdfe)
        mu <- pmax(exp(eta), eps100)
        is_sh <- TRUE; ok <- 0L; num_sh <- num_sh + 1L
      } else {
        ok <- 0L; num_sh <- 0L
      }
    }

    if (verbose > -1L) {
      cat(sprintf("Iter %3d:  dev = %-11.5e", it, deviance))
      if (it > 1L) cat(sprintf("  eps = %-9.4e", eps_val))
      if (K > 0L && is.finite(beta_change)) cat(sprintf("  db = %-9.4e", beta_change))
      cat(sprintf("  tol = %5.0e", inner_tol))
      if (is_sh) cat("  H")
      if (ok > 0L) cat("  O")
      cat("\n")
    }

    if (is_sh) { deviance <- old_deviance; next }
    if (converged) break

    # Adaptive inner tolerance
    if (it > 1L && eps_val < inner_tol) {
      inner_tol <- max(min(0.1 * inner_tol, alt_tol), target_inner_tol)
      alt_tol <- 10^(-ceiling(log10(1 / max(0.1 * eps_val, .Machine$double.eps))))
    }
  }

  if (!converged) {
    stop(sprintf(paste0(
      "IRLS failed to converge in %d iterations (eps = %.4e); ",
      "coefficients are numerically meaningless and will not be returned. ",
      "Try increasing maxiter, loosening tol, standardize = TRUE, or a ",
      "simpler FE structure."), maxiter, eps_val), call. = FALSE)
  }
  if (verbose > -1L) {
    cat(sprintf("Converged in %d iterations (tol = %9.4e)\n", iter, tol))
  }

  # ==================================================================
  # Final beta and VCE
  # ==================================================================
  irls_w <- pmax(w_user * mu, 1e-20)
  offv <- offset_orig             # restore original offset scale
  swf <- sum(irls_w)
  mean_eta_no_off <- sum(irls_w * (eta - offv)) / swf
  mean_X <- colSums(irls_w * X) / swf
  b_cons <- mean_eta_no_off - drop(crossprod(mean_X, b))

  # Recompute Xhat for the VCE
  ZwZf <- crossprod(Z_dm, irls_w * Z_dm)
  Pif <- solve(ZwZf, crossprod(Z_dm, irls_w * X_dm))
  Xhat <- Z_dm %*% Pif

  # Effective N and meat weight (fweight convention from ppmlhdfe)
  N_eff <- if (weight_type == "fweight") sum(w_user) else N
  w_meat <- if (weight_type == "fweight") sqrt(w_user) * mu else irls_w

  bread <- tryCatch(solve(crossprod(Xhat, irls_w * X_dm)),
                    error = function(e) NULL)
  actual_rank <- if (is.null(bread)) {
    qr(crossprod(Xhat, irls_w * X_dm))$rank
  } else K

  G_counts <- integer(0)
  vcetype <- "robust"
  if (has_cluster) {
    vcetype <- "cluster"
    vc <- vce_cluster(Xhat, X_dm, irls_w, w_meat, resid, clust_df, N_eff)
    V_slope <- vc$V
    G_counts <- vc$G_counts
  } else {
    V_slope <- vce_robust(Xhat, X_dm, irls_w, w_meat, resid, N_eff)
  }

  # Back-transform standardized coefficients and VCE
  if (standardize && K > 0L) {
    b <- b / stdev_x
    V_slope <- V_slope / tcrossprod(stdev_x)
    X <- sweep(X, 2L, stdev_x, "*")
    mean_X <- colSums(irls_w * X) / swf
    b_cons <- mean_eta_no_off - drop(crossprod(mean_X, b))
  }
  b_full <- c(b, b_cons)

  # Expand V to include _cons (zero VCE: partialled out, ivreg2 convention)
  K_total <- K + 1L
  V <- matrix(0, K_total, K_total)
  V[seq_len(K), seq_len(K)] <- V_slope

  # ---- Log pseudo-likelihood ----------------------------------------------
  ll <- sum(w_user * (y * eta - mu - lgamma(y + 1)))
  if (has_offset && sum(abs(offv)) > 0) {
    ll0c <- log(sum(w_user * y) / sum(w_user * exp(offv)))
    mu0 <- exp(ll0c + offv)
    ll_0 <- sum(w_user * (y * (ll0c + offv) - mu0 - lgamma(y + 1)))
  } else {
    ll_0 <- sum(w_user * (y * log(mean_y) - mean_y - lgamma(y + 1)))
  }

  # ---- d values: FE sum ----------------------------------------------------
  d_vals <- eta - offv
  if (K > 0L) d_vals <- d_vals - drop(X %*% b)
  d_vals <- d_vals - b_cons
  d_vals <- d_vals - sum(irls_w * d_vals) / sum(irls_w)   # weighted-mean-0
  d_vals <- unname(d_vals)

  # ---- Degrees of freedom ---------------------------------------------------
  df_a <- compute_df_a(fe_ids,
                       cluster_ids = if (has_cluster) as.list(clust_df) else NULL)
  if (has_cluster) {
    df_resid <- min(G_counts) - 1L
    N_clust <- min(G_counts)
  } else {
    df_resid <- N_eff - K - df_a
    N_clust <- NA_integer_
  }

  # Wald chi2 on the slope block
  chi2 <- tryCatch(
    drop(t(b) %*% solve(V_slope, b)),
    error = function(e) NA_real_)

  names(b_full) <- c(colnames_x, "_cons")
  dimnames(V) <- list(names(b_full), names(b_full))

  out <- list(
    coefficients = b_full,
    vcov = V,
    coefnames = names(b_full),
    exog_names = if (K_exog) colnames_x[seq_len(K_exog)] else character(0),
    endog_names = if (K_endo) colnames_x[K_exog + seq_len(K_endo)] else character(0),
    inst_names = inst_names,
    N = N_eff,
    N_rows = N,
    N_full = N_input,
    converged = converged,
    iterations = iter,
    deviance = deviance,
    ll = ll,
    ll_0 = ll_0,
    r2_p = 1 - ll / ll_0,
    chi2 = chi2,
    df_m = actual_rank,
    rank = actual_rank,
    df = df_resid,
    df_a = df_a,
    vcetype = vcetype,
    N_clust = N_clust,
    G_counts = G_counts,
    weight_type = if (has_weight) weight_type else "",
    num_singletons = num_singletons,
    num_sep_fe = num_sep_fe,
    num_sep_advanced = num_sep_advanced,
    num_sep_mu = num_sep_mu,
    num_separated = num_sep_fe + num_sep_advanced + num_sep_mu,
    n_dropped_exog = n_dropped_exog,
    n_dropped_endog = n_dropped_endog,
    n_dropped_inst = n_dropped_inst,
    fitted.values = mu,
    eta = eta,
    d = d_vals,
    residuals = y - mu,
    x = if (K > 0L) `colnames<-`(X, colnames_x) else matrix(0, N, 0L),
    y = y,
    offset = if (has_offset) offv else NULL,
    esample = esample,
    sep_mask = sep_mask,
    tol = tol,
    call = NULL,
    fml = NULL
  )
  class(out) <- "ivppmlhdfe"
  out
}