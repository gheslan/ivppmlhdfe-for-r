# spj.R -- split-panel jackknife bias correction with bootstrap SE
#
# Implements the SPJ formulas of Kwon, Larch, Yoon & Yotov (2026),
# Section 2.3 / Appendix B, matching the template do-files distributed with
# the Stata package (MC_SPJ_BTS_ClassA/B/C.do):
#
#   Class A:  beta_SPJ = 3 b_full - mean(b_T1, b_T2) - mean(b_N1, b_N2)
#   Class B:  beta_SPJ = 2 b_full - mean(b_aa, b_ab, b_ba, b_bb)
#   Class C:  beta_SPJ = 4 b_full - 2 mean(country) - 2 mean(time)
#                        + mean(8 cells)
#   Interacted single-family FE: beta_SPJ = 2 b_full - mean(b_N1, b_N2)
#
# Bootstrap resampling follows the paper: individuals (Class A) or directed
# pairs (Classes B and C) are resampled with replacement, retaining all T
# periods per sampled unit.

#' Split-panel jackknife bias correction for IV-PPML-HDFE
#'
#' Computes the SPJ-corrected point estimate of the parameter(s) of interest
#' for the three fixed-effect structures studied in Kwon, Larch, Yoon and
#' Yotov (2026), optionally with pair/individual-bootstrap standard errors.
#'
#' @param fml,data,... passed to [ivppmlhdfe()] for every (sub-)estimation.
#' @param class one of `"A"` (individual + time FE; splits along both the
#'   time and the cross-sectional dimension), `"B"` (two-way gravity FE;
#'   4-subpanel country split), `"C"` (three-way gravity FE; 8-panel
#'   country-by-time split), or `"interacted"` (single interacted FE family;
#'   cross-sectional split only, as in the Chen et al. (2025) replication).
#' @param id (Class A / interacted) column name of the individual id.
#' @param time (Classes A and C) column name of the time variable.
#' @param exp,imp (Classes B and C) column names of the origin/destination
#'   country ids.
#' @param coef_name name(s) of the coefficient(s) to bias-correct; default
#'   is the (first) endogenous regressor.
#' @param B number of bootstrap replications (0 = no bootstrap; the paper
#'   uses `B = 1000`).
#' @param cluster_boot column name of the bootstrap resampling unit;
#'   defaults to `id` for Classes A/interacted and the directed pair
#'   (`exp` x `imp`) for Classes B/C.
#' @param seed optional RNG seed.
#' @param verbose print progress.
#'
#' @return A list of class `"ivppml_spj"`:
#'   \item{b_full}{uncorrected full-sample estimate.}
#'   \item{b_spj}{SPJ-corrected estimate.}
#'   \item{subpanel}{named vector/list of sub-panel estimates.}
#'   \item{boot}{matrix of bootstrap SPJ estimates (successful draws).}
#'   \item{se_boot}{bootstrap standard deviation.}
#'   \item{se_ci_implied}{the paper's preferred aggregator: width of the
#'     95 percent percentile interval divided by 2 x 1.96 (robust to a small
#'     fraction of divergent draws).}
#'   \item{ci}{95 percent percentile confidence interval.}
#'   \item{n_failed}{number of bootstrap draws on which any sub-estimation
#'     failed (these are skipped, matching the `capture` logic of the
#'     template do-files).}
#'
#' @examples
#' \donttest{
#' data(ivppmlhdfe_ClassA)
#' s <- ivppml_spj(y ~ x2 | id + year | x1 ~ z, data = ivppmlhdfe_ClassA,
#'                 class = "A", id = "id", time = "year", B = 0)
#' s$b_spj
#' }
#' @export
ivppml_spj <- function(fml, data, class = c("A", "B", "C", "interacted"),
                       id = NULL, time = NULL, exp = NULL, imp = NULL,
                       coef_name = NULL, B = 0L,
                       cluster_boot = NULL, seed = NULL,
                       verbose = FALSE, ...) {
  class <- match.arg(class)
  if (!is.null(seed)) set.seed(seed)
  data <- as.data.frame(data)

  if (class %in% c("A", "interacted") && is.null(id))
    stop("class \"", class, "\" requires `id`")
  if (class %in% c("A", "C") && is.null(time))
    stop("class \"", class, "\" requires `time`")
  if (class %in% c("B", "C") && (is.null(exp) || is.null(imp)))
    stop("class \"", class, "\" requires `exp` and `imp`")

  est <- function(d) {
    fit <- ivppmlhdfe(fml, data = d, ...)
    bb <- coef(fit)
    nm <- if (is.null(coef_name)) {
      if (length(fit$endog_names)) fit$endog_names else names(bb)[1L]
    } else coef_name
    bb[nm]
  }

  spj_once <- function(d) {
    b_full <- est(d)
    sub <- list()
    if (class == "A") {
      tmid <- floor((min(d[[time]]) + max(d[[time]])) / 2)
      bT1 <- est(d[d[[time]] <= tmid, ]); bT2 <- est(d[d[[time]] > tmid, ])
      ids <- unique(d[[id]])
      half <- stats::setNames(stats::runif(length(ids)) < 0.5, ids)
      h <- half[as.character(d[[id]])]
      bN1 <- est(d[h, ]); bN2 <- est(d[!h, ])
      b_spj <- 3 * b_full - (bT1 + bT2) / 2 - (bN1 + bN2) / 2
      sub <- list(T1 = bT1, T2 = bT2, N1 = bN1, N2 = bN2)
    } else if (class == "interacted") {
      ids <- unique(d[[id]])
      half <- stats::setNames(stats::runif(length(ids)) < 0.5, ids)
      h <- half[as.character(d[[id]])]
      bN1 <- est(d[h, ]); bN2 <- est(d[!h, ])
      b_spj <- 2 * b_full - (bN1 + bN2) / 2
      sub <- list(N1 = bN1, N2 = bN2)
    } else {
      countries <- unique(c(d[[exp]], d[[imp]]))
      half <- stats::setNames(stats::runif(length(countries)) < 0.5, countries)
      he <- half[as.character(d[[exp]])]
      hi <- half[as.character(d[[imp]])]
      b_cells <- numeric(0)
      for (ee in c(FALSE, TRUE)) for (ii in c(FALSE, TRUE)) {
        b_cells <- c(b_cells, est(d[he == ee & hi == ii, ]))
      }
      if (class == "B") {
        b_spj <- 2 * b_full - mean(b_cells)
        sub <- list(country = b_cells)
      } else {
        tmid <- floor((min(d[[time]]) + max(d[[time]])) / 2)
        ht <- d[[time]] <= tmid
        b_t <- c(est(d[ht, ]), est(d[!ht, ]))
        b8 <- numeric(0)
        for (ee in c(FALSE, TRUE)) for (ii in c(FALSE, TRUE))
          for (tt in c(TRUE, FALSE)) {
            b8 <- c(b8, est(d[he == ee & hi == ii & (ht == tt), ]))
          }
        b_spj <- 4 * b_full - 2 * mean(b_cells) - 2 * mean(b_t) + mean(b8)
        sub <- list(country = b_cells, time = b_t, cell8 = b8)
      }
    }
    list(b_full = b_full, b_spj = b_spj, subpanel = sub)
  }

  main <- spj_once(data)

  boot_mat <- NULL; n_failed <- 0L
  se_boot <- NA_real_; se_ci <- NA_real_; ci <- c(NA_real_, NA_real_)
  if (B > 0L) {
    # default bootstrap resampling unit
    if (is.null(cluster_boot)) {
      boot_id <- if (class %in% c("A", "interacted")) data[[id]]
                 else interact_id(list(data[[exp]], data[[imp]]))
    } else {
      boot_id <- data[[cluster_boot]]
    }
    groups <- split(seq_len(nrow(data)), boot_id)
    G <- length(groups)
    draws <- vector("list", B)
    for (bb in seq_len(B)) {
      pick <- sample.int(G, G, replace = TRUE)
      rows <- unlist(groups[pick], use.names = FALSE)
      dB <- data[rows, , drop = FALSE]
      # relabel resampled clusters (bsample idcluster): each draw is a new id
      reps <- vapply(groups[pick], length, integer(1))
      newid <- rep(seq_len(G), reps)
      if (class %in% c("A", "interacted")) {
        dB[[id]] <- newid
      } else {
        # new directed-pair id; exp/imp labels are retained for the splits
        dB[["..pair_boot"]] <- newid
      }
      r <- tryCatch(spj_once(dB)$b_spj, error = function(e) NULL)
      if (is.null(r) || any(!is.finite(r))) { n_failed <- n_failed + 1L; next }
      draws[[bb]] <- r
      if (verbose && bb %% 50L == 0L) {
        cat(sprintf("bootstrap %d/%d (failed: %d)\n", bb, B, n_failed))
      }
    }
    ok <- !vapply(draws, is.null, logical(1))
    if (any(ok)) {
      boot_mat <- do.call(rbind, draws[ok])
      se_boot <- apply(boot_mat, 2L, stats::sd)
      qs <- apply(boot_mat, 2L, stats::quantile, probs = c(0.025, 0.975))
      ci <- t(qs)
      se_ci <- (qs[2L, ] - qs[1L, ]) / (2 * stats::qnorm(0.975))
    }
  }

  out <- list(b_full = main$b_full, b_spj = main$b_spj,
              subpanel = main$subpanel, class = class,
              boot = boot_mat, se_boot = se_boot,
              se_ci_implied = se_ci, ci = ci, n_failed = n_failed, B = B)
  class(out) <- "ivppml_spj"
  out
}

#' @export
print.ivppml_spj <- function(x, ...) {
  cat("Split-panel jackknife (Class ", x$class, ")\n", sep = "")
  cat("  Full-sample estimate: ", format(x$b_full), "\n")
  cat("  SPJ estimate:         ", format(x$b_spj), "\n")
  if (!is.null(x$boot)) {
    cat("  Bootstrap SE (sd):    ", format(x$se_boot), "\n")
    cat("  CI-implied SE:        ", format(x$se_ci_implied), "\n")
    cat("  95% percentile CI:    [", format(x$ci[1L]), ", ",
        format(x$ci[2L]), "]\n", sep = "")
    cat("  Successful draws:     ", nrow(x$boot), "/", x$B,
        " (", x$n_failed, " failed)\n", sep = "")
  }
  invisible(x)
}