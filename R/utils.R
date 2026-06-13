# utils.R -- internal utilities for ivppmlhdfe
#
# These helpers mirror low-level building blocks of the Stata/Mata
# implementation (reghdfe's fixed_effects() machinery, group-id creation,
# and the weighted within-transformation used at every IRLS iteration).

#' @importFrom stats ave complete.cases pchisq pnorm qnorm quantile sd setNames
NULL

# ---------------------------------------------------------------------------
# Group ids: integer codes 1..G for a single factor-like vector, or the
# interaction of several vectors (Stata's `egen group(...)` / `i.a#i.b`).
# ---------------------------------------------------------------------------
group_id <- function(...) {
  args <- list(...)
  if (length(args) == 1L) {
    x <- args[[1L]]
    if (is.factor(x)) return(as.integer(x))
    return(match(x, unique(x)))
  }
  key <- do.call(paste, c(args, sep = "\r"))
  match(key, unique(key))
}

# Interaction id from a list of vectors (used for `a^b` FE terms and
# interaction clusters). Equivalent to Mata's ivppmlhdfe_interact_id chains.
interact_id <- function(vars) {
  do.call(group_id, vars)
}

# ---------------------------------------------------------------------------
# Weighted multi-way demeaning (the inner HDFE solve).
#
# Wraps fixest::demean(), which implements the same alternating-projection
# (method of alternating projections with acceleration) algorithm family as
# reghdfe's partial_out. `tol` maps to the adaptive inner tolerance of the
# Stata implementation (HDFE.tolerance).
# ---------------------------------------------------------------------------
hdfe_demean <- function(M, fe_df, weights, tol) {
  if (is.null(fe_df) || ncol(fe_df) == 0L) {
    # noabsorb: demean by the constant only (weighted grand mean)
    sw <- sum(weights)
    mns <- colSums(M * weights) / sw
    return(sweep(M, 2L, mns, "-"))
  }
  out <- fixest::demean(M, fe_df, weights = weights,
                        tol = tol, iter = 100000L,
                        notes = FALSE, im_confident = TRUE)
  as.matrix(out)
}

# ---------------------------------------------------------------------------
# edittozerotol: round |x| < tol to exact zero (matches Mata's
# _edittozerotol applied before/after partialling out).
# ---------------------------------------------------------------------------
edit_to_zero <- function(M, tol) {
  M[abs(M) < tol] <- 0
  M
}

# ---------------------------------------------------------------------------
# Parse a cluster specification into a data.frame of id vectors.
# Accepts: NULL, a one-sided formula (~c1 + c2, with `^` interactions),
# a character vector of column names (names may contain "#" or "^" for
# interactions, matching Stata's `vce(cluster exp#imp)`), or a vector /
# data.frame of ids supplied directly.
# ---------------------------------------------------------------------------
parse_cluster <- function(cluster, data, n) {
  if (is.null(cluster)) return(NULL)
  if (inherits(cluster, "formula")) {
    terms_chr <- attr(stats::terms(cluster), "term.labels")
    cols <- lapply(terms_chr, function(tt) eval_fe_term(tt, data))
    out <- as.data.frame(cols, optional = TRUE)
    names(out) <- terms_chr
    return(out)
  }
  if (is.character(cluster)) {
    cols <- lapply(cluster, function(tt) eval_fe_term(tt, data))
    out <- as.data.frame(cols, optional = TRUE)
    names(out) <- cluster
    return(out)
  }
  if (is.data.frame(cluster)) {
    stopifnot(nrow(cluster) == n)
    return(cluster)
  }
  if (is.atomic(cluster)) {
    stopifnot(length(cluster) == n)
    return(data.frame(cluster1 = cluster))
  }
  if (is.list(cluster)) {
    out <- as.data.frame(cluster, optional = TRUE)
    stopifnot(nrow(out) == n)
    return(out)
  }
  stop("invalid `cluster` specification")
}

# Evaluate one FE / cluster term such as "id", "exp^year", or "exp#year"
# against `data`, returning an id vector.
eval_fe_term <- function(term, data) {
  term <- gsub("`", "", term, fixed = TRUE)
  term <- gsub("#", "^", term, fixed = TRUE)
  parts <- strsplit(term, "^", fixed = TRUE)[[1L]]
  parts <- trimws(parts)
  vars <- lapply(parts, function(p) {
    v <- tryCatch(eval(parse(text = p), envir = data, enclos = parent.frame()),
                  error = function(e) NULL)
    if (is.null(v)) stop("variable `", p, "` not found in `data`", call. = FALSE)
    v
  })
  if (length(vars) == 1L) vars[[1L]] else interact_id(vars)
}

# ---------------------------------------------------------------------------
# Absorbed-FE degrees of freedom (reghdfe's estimate_dof, default dof(all)):
#   df_a = G1 + (G2 - M2) + sum_{k >= 3} (G_k - 1)
# where M2 is the number of connected components ("mobility groups") of the
# bipartite graph linking the first two FEs. FEs nested within a cluster
# variable contribute 0 (reghdfe convention). Only affects e(df) / rmse,
# never coefficients or VCE.
# ---------------------------------------------------------------------------
n_connected_components <- function(id1, id2) {
  # union-find over levels of id1 (1..G1) and id2 (G1+1..G1+G2)
  g1 <- max(id1); g2 <- max(id2)
  parent <- seq_len(g1 + g2)
  find <- function(i) {
    while (parent[i] != i) {
      parent[i] <<- parent[parent[i]]
      i <- parent[i]
    }
    i
  }
  # iterate over unique (id1, id2) pairs only
  key <- id1 * (g2 + 1) + id2
  keep <- !duplicated(key)
  a <- id1[keep]; b <- id2[keep] + g1
  for (k in seq_along(a)) {
    ra <- find(a[k]); rb <- find(b[k])
    if (ra != rb) parent[rb] <- ra
  }
  roots <- vapply(seq_len(g1 + g2), find, integer(1))
  length(unique(roots[unique(c(id1, id2 + g1))]))
}

compute_df_a <- function(fe_ids, cluster_ids = NULL) {
  n_fe <- length(fe_ids)
  if (n_fe == 0L) return(0L)
  # nested-in-cluster check: FE k is nested in cluster c if each FE level
  # maps to a single cluster level
  nested <- rep(FALSE, n_fe)
  if (!is.null(cluster_ids) && length(cluster_ids)) {
    for (k in seq_len(n_fe)) {
      for (cid in cluster_ids) {
        tab <- rowsum(rep(1L, length(cid)),
                      group = paste(fe_ids[[k]], group_id(cid), sep = "\r"))
        n_pairs <- nrow(tab)
        if (n_pairs == length(unique(fe_ids[[k]]))) { nested[k] <- TRUE; break }
      }
    }
  }
  levs <- vapply(fe_ids, function(v) length(unique(v)), integer(1))
  df <- integer(n_fe)
  df[1L] <- levs[1L]
  if (n_fe >= 2L) {
    M2 <- n_connected_components(group_id(fe_ids[[1L]]), group_id(fe_ids[[2L]]))
    df[2L] <- levs[2L] - M2
  }
  if (n_fe >= 3L) for (k in 3:n_fe) df[k] <- levs[k] - 1L
  df[nested] <- 0L
  sum(df)
}

# ---------------------------------------------------------------------------
# Weighted column standard deviations (population variance, matching Mata's
# quadvariance + the Julia implementation's dot(w, (x - mx)^2)/sum(w)).
# ---------------------------------------------------------------------------
wcolsd <- function(M, w) {
  sw <- sum(w)
  apply(M, 2L, function(x) {
    mx <- sum(w * x) / sw
    sqrt(max(sum(w * (x - mx)^2) / sw, 0))
  })
}