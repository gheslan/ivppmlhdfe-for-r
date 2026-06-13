# collinearity.R -- two-stage collinearity removal and VCE estimators
#
# Mirrors ivppmlhdfe.ado v0.9.4 (Mata) and the Julia backend v0.9.4.

# ---------------------------------------------------------------------------
# _select_not_collinear: rank detection on FE-partialled columns with an
# absorption test against the raw (pre-partial) columns.
#
# Two-part test per column j (verbatim from the Julia reference):
#   (a) Absorption: weighted norm of the partialled column must exceed
#       rtol * max(raw weighted norms). Catches columns fully absorbed by
#       the FE structure.
#   (b) Within-block rank: pivoted QR on the survivors with relative
#       threshold rtol * max|R_ii|.
#
# Returns sorted integer indices of surviving columns.
# ---------------------------------------------------------------------------
select_not_collinear <- function(A_partial, A_raw, w, rtol = 1e-9) {
  K <- ncol(A_partial)
  if (K == 0L) return(integer(0))
  sw <- sqrt(pmax(w, 0))
  raw_norms <- sqrt(colSums((sw * A_raw)^2))
  max_raw <- max(raw_norms, 0)
  if (max_raw == 0) return(integer(0))
  abs_thresh <- rtol * max_raw

  partial_norms <- sqrt(colSums((sw * A_partial)^2))
  nonzero_idx <- which(partial_norms > abs_thresh)
  if (length(nonzero_idx) <= 1L) return(nonzero_idx)

  WA <- sw * A_partial[, nonzero_idx, drop = FALSE]
  qrf <- qr(WA, LAPACK = TRUE)         # pivoted (column-norm) QR
  diagR <- abs(diag(qr.R(qrf)))
  thresh <- rtol * max(diagR)
  kept_pivoted <- which(diagR > thresh)
  kept_in_sub <- qrf$pivot[kept_pivoted]
  sort(nonzero_idx[kept_in_sub])
}

# ---------------------------------------------------------------------------
# Robust sandwich VCE (ppmlhdfe convention: N/(N-1) multiplier).
#   w_bread: bread weighting (irls_w = w_user * mu in all cases)
#   w_meat : per-row score weighting:
#            pw/aw/unweighted -> irls_w        (meat scales as w_user^2)
#            fweight          -> sqrt(w) * mu  (meat scales as w_user)
# ---------------------------------------------------------------------------
vce_robust <- function(Xhat, X_dm, w_bread, w_meat, resid, N_eff) {
  bread <- solve(crossprod(Xhat, w_bread * X_dm))
  s <- Xhat * (w_meat * resid)
  meat <- crossprod(s)
  (N_eff / (N_eff - 1)) * bread %*% meat %*% bread
}

# Cluster meat: sum of outer products of within-cluster score sums.
clust_meat <- function(scores, clust_id) {
  K <- ncol(scores)
  agg <- rowsum(scores, group = clust_id, reorder = FALSE)
  list(meat = crossprod(agg), G = nrow(agg))
}

# ---------------------------------------------------------------------------
# Cluster-robust VCE. Single cluster uses the Arellano G/(G-1) correction;
# multi-way uses the Cameron-Gelbach-Miller inclusion-exclusion formula with
# the reghdfe PSD fix (shift eigenvalues if min eigenvalue < 0).
# Returns list(V, G_counts).
# ---------------------------------------------------------------------------
vce_cluster <- function(Xhat, X_dm, w_bread, w_meat, resid, clust_df, N_eff) {
  K <- ncol(X_dm)
  bread <- solve(crossprod(Xhat, w_bread * X_dm))
  scores <- Xhat * (w_meat * resid)
  n_clust <- ncol(clust_df)

  if (n_clust == 1L) {
    cm <- clust_meat(scores, clust_df[[1L]])
    if (cm$G <= 1L) {
      warning("only ", cm$G, " cluster(s); cluster VCE not defined. Using robust.")
      V <- vce_robust(Xhat, X_dm, w_bread, w_meat, resid, N_eff)
    } else {
      V <- (cm$G / (cm$G - 1)) * bread %*% cm$meat %*% bread
    }
    return(list(V = V, G_counts = cm$G))
  }

  V <- matrix(0, K, K)
  G_counts <- integer(n_clust)
  for (mask in seq_len(2^n_clust - 1L)) {
    in_set <- which(bitwAnd(mask, bitwShiftL(1L, seq_len(n_clust) - 1L)) != 0L)
    cid <- clust_df[[in_set[1L]]]
    if (length(in_set) > 1L) {
      for (j in in_set[-1L]) cid <- interact_id(list(cid, clust_df[[j]]))
    }
    cm <- clust_meat(scores, cid)
    sgn <- if (length(in_set) %% 2L == 1L) 1 else -1
    if (cm$G > 1L) {
      V <- V + sgn * (cm$G / (cm$G - 1)) * bread %*% cm$meat %*% bread
    }
    if (length(in_set) == 1L) G_counts[in_set] <- cm$G
  }
  V <- 0.5 * (V + t(V))                 # symmetrize
  ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < 0) V <- V - min(ev) * diag(K)   # PSD fix (reghdfe)
  list(V = V, G_counts = G_counts)
}