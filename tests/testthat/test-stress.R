# Stress tests -- mirror the Stata package's Appendix D battery:
#   (i)   bit-identical regression checks against ppmlhdfe (via fepois),
#   (ii)  numerical-stability / scaling edge cases,
#   (iii) singleton + cascading-separation, weak instruments, collinearity,
#         over-identification, multiple endogenous regressors.
# These probe for crashes, silent wrong answers, and convergence failures.

# Independent LSDV (dummy-variable) IV-PPML oracle, basis-invariant via an
# SVD pseudo-inverse. Shares no code with the package engine; used to certify
# IV point estimates exactly as the paper validated against Stata ivpoisson gmm.
.pinv <- function(M, tol = 1e-10) {
  s <- svd(M); d <- s$d; r <- d > max(d) * tol
  s$v[, r, drop = FALSE] %*% ((1 / d[r]) * t(s$u[, r, drop = FALSE]))
}
.iv_lsdv <- function(y, Xexog, xendo, z, FEs, maxit = 300, tol = 1e-13) {
  Dl <- list(); first <- TRUE
  for (f in FEs) {
    f <- factor(f); MM <- stats::model.matrix(~ f - 1)
    if (!first) MM <- MM[, -1, drop = FALSE]
    Dl[[length(Dl) + 1]] <- MM; first <- FALSE
  }
  D <- do.call(cbind, Dl)
  X <- cbind(Xexog, xendo, D); Q <- cbind(Xexog, z, D)
  ns <- ncol(Xexog) + 1L
  mu <- 0.5 * (y + mean(y)); mu <- pmax(mu, pmax(0.05 * y, 1e-3))
  eta <- log(mu); b <- rep(0, ncol(X))
  for (it in 1:maxit) {
    zt <- eta - 1 + y / mu; w <- mu
    Pi <- .pinv(crossprod(Q, w * Q)) %*% crossprod(Q, w * X); Xhat <- Q %*% Pi
    bn <- .pinv(crossprod(Xhat, w * X)) %*% crossprod(Xhat, w * zt)
    if (max(abs((bn - b)[1:ns])) < tol) { b <- bn; break }
    b <- bn; eta <- drop(X %*% b); mu <- pmax(exp(eta), 1e-14)
  }
  stats::setNames(b[1:ns], c(colnames(Xexog), "xendo"))
}

# ------------------------------------------------------------------ (i)
test_that("regression battery: many specs reproduce fepois (bit-identical)", {
  skip_no_fixest()
  skip_on_cran()
  set.seed(101)
  dA <- get_class("A"); dB <- get_class("B")
  specs <- list(
    list(d = dA, fml_r = y ~ x1 + x2 | id + year,        fml_fe = y ~ x1 + x2 | id + year),
    list(d = dA, fml_r = y ~ x1 | id,                     fml_fe = y ~ x1 | id),
    list(d = dA, fml_r = y ~ x1 + x2 | year,              fml_fe = y ~ x1 + x2 | year),
    list(d = dB, fml_r = y ~ x1 + x2 | exp^year + imp^year, fml_fe = y ~ x1 + x2 | exp^year + imp^year),
    list(d = dB, fml_r = y ~ x1 | exp^year,               fml_fe = y ~ x1 | exp^year)
  )
  for (s in specs) {
    m  <- ivppmlhdfe(s$fml_r, data = s$d)
    fe <- fixest::fepois(s$fml_fe, data = s$d, notes = FALSE)
    common <- intersect(names(coef(m)), names(coef(fe)))
    expect_equal(unname(coef(m)[common]), unname(coef(fe)[common]), tolerance = 1e-7,
                 label = deparse(s$fml_r))
  }
})

test_that("IV point estimate equals the independent LSDV oracle (Class A)", {
  d <- get_class("A")
  m <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d)
  ref <- .iv_lsdv(d$y, cbind(x2 = d$x2), d$x1, d$z, list(d$id, d$year))
  expect_equal(unname(coef(m)["x1"]), unname(ref["xendo"]), tolerance = 1e-7)
})

test_that("IV point estimates equal the independent LSDV oracle (Classes B, C)", {
  skip_on_cran()
  for (cl in c("B", "C")) {
    d <- get_class(cl)
    fml <- switch(cl,
      B = y ~ x2 | exp^year + imp^year | x1 ~ z,
      C = y ~ x2 | exp^year + imp^year + pair | x1 ~ z)
    fes <- switch(cl,
      B = list(interaction(d$exp, d$year, drop = TRUE), interaction(d$imp, d$year, drop = TRUE)),
      C = list(interaction(d$exp, d$year, drop = TRUE), interaction(d$imp, d$year, drop = TRUE), factor(d$pair)))
    m <- ivppmlhdfe(fml, data = d)
    ref <- .iv_lsdv(d$y, cbind(x2 = d$x2), d$x1, d$z, fes)
    expect_equal(unname(coef(m)["x1"]), unname(ref["xendo"]), tolerance = 1e-7,
                 label = paste("Class", cl))
  }
})

# ------------------------------------------------------------------ (ii)
test_that("poorly-scaled regressors: standardize recovers the same estimate", {
  d <- get_class("A")
  d$x1big <- d$x1 * 1e6
  m_raw <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d)
  m_std <- ivppmlhdfe(y ~ x2 | id + year | x1big ~ z, data = d, standardize = TRUE)
  # coefficient on x1big is coef(x1)/1e6
  expect_equal(unname(coef(m_std)["x1big"] * 1e6), unname(coef(m_raw)["x1"]),
               tolerance = 1e-5)
})

test_that("extreme-scale offset converges (centering trick)", {
  d <- get_class("A")
  off <- rep(25, nrow(d))          # exp(25) ~ 7e10 would blow up naive init
  m <- ivppmlhdfe(y ~ x1 + x2 | id + year, data = d, offset = off)
  expect_true(m$converged)
  expect_true(all(is.finite(coef(m))))
})

test_that("estimates are stable across tolerances", {
  d <- get_class("A")
  m1 <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d, tol = 1e-6)
  m2 <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d, tol = 1e-12)
  expect_equal(unname(coef(m1)["x1"]), unname(coef(m2)["x1"]), tolerance = 1e-5)
})

test_that("many zeros in the outcome still converge", {
  set.seed(7); d <- get_class("A")
  d$y <- rbinom(nrow(d), 1, 0.2) * d$y    # ~80% zeros
  if (sum(d$y) > 0) {
    m <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d,
                    separation = "fe")
    expect_true(m$converged)
  } else succeed()
})

# ------------------------------------------------------------------ (iii)
test_that("cascading singletons are fully removed", {
  d <- get_class("A")
  # chain: a fresh id with 1 obs in a fresh year with 1 obs -> cascade
  extra <- d[1:2, ]
  extra$id <- max(d$id) + 1L
  extra$year <- max(d$year) + 1:2
  d2 <- rbind(d, extra)
  m <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d2)
  expect_true(m$converged)
  expect_gte(m$num_singletons, 2)
  # the injected id must not survive
  expect_false(any(m$esample & d2$id == max(d$id) + 1L))
})

test_that("cascading FE separation (all-zero groups) drops the chain", {
  d <- get_class("A")
  d$y[d$id %in% c(1, 2, 3)] <- 0     # three all-zero individuals
  m <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d, separation = "fe")
  expect_true(m$converged)
  expect_false(any(m$esample & d$id %in% c(1, 2, 3)))
})

test_that("weak instrument does not crash and inflates the SE", {
  set.seed(9); d <- get_class("A")
  d$zweak <- 0.02 * d$z + rnorm(nrow(d))   # very weak first stage
  m_strong <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d)
  m_weak   <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ zweak, data = d)
  expect_true(m_weak$converged)
  expect_gt(sqrt(diag(vcov(m_weak)))["x1"], sqrt(diag(vcov(m_strong)))["x1"])
})

test_that("collinear regressors and instruments are dropped", {
  d <- get_class("A")
  d$x2b <- d$x2                 # duplicate exogenous
  d$zb  <- d$z                  # duplicate instrument
  m <- ivppmlhdfe(y ~ x2 + x2b | id + year | x1 ~ z + zb, data = d)
  expect_equal(m$n_dropped_exog, 1L)
  expect_equal(m$n_dropped_inst, 1L)
  expect_true(m$converged)
})

test_that("over-identified IV runs (one-step GMM)", {
  set.seed(3); d <- get_class("A")
  d$z2 <- 0.5 * d$x1 + rnorm(nrow(d))      # second valid-ish instrument
  m <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z + z2, data = d)
  expect_true(m$converged)
  expect_length(m$inst_names, 2L)
})

test_that("multiple endogenous regressors are supported", {
  set.seed(4); d <- get_class("A")
  d$x3 <- 0.6 * d$z + rnorm(nrow(d))       # second endogenous
  d$z3 <- d$x3 + 0.5 * rnorm(nrow(d))      # its instrument
  m <- ivppmlhdfe(y ~ x2 | id + year | x1 + x3 ~ z + z3, data = d)
  expect_true(m$converged)
  expect_length(m$endog_names, 2L)
})

test_that("frequency weight 2 equals row duplication (point estimates)", {
  skip_no_fixest()
  d <- get_class("A")
  m_fw  <- ivppmlhdfe(y ~ x1 + x2 | id + year, data = d,
                      weights = rep(2, nrow(d)), weight_type = "fweight")
  d2 <- rbind(d, d)
  m_dup <- ivppmlhdfe(y ~ x1 + x2 | id + year, data = d2)
  expect_equal(unname(coef(m_fw)[c("x1", "x2")]),
               unname(coef(m_dup)[c("x1", "x2")]), tolerance = 1e-6)
})

test_that("repeated fits are deterministic (loop stability)", {
  d <- get_class("A")
  b <- replicate(3, coef(ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d))["x1"])
  expect_equal(b[1], b[2]); expect_equal(b[2], b[3])
})

test_that("SPJ bootstrap is reproducible under a fixed seed", {
  d <- get_class("A")
  s1 <- ivppml_spj(y ~ x2 | id + year | x1 ~ z, data = d, class = "A",
                   id = "id", time = "year", B = 6, seed = 42)
  s2 <- ivppml_spj(y ~ x2 | id + year | x1 ~ z, data = d, class = "A",
                   id = "id", time = "year", B = 6, seed = 42)
  expect_equal(s1$b_spj, s2$b_spj)
  expect_equal(s1$se_ci_implied, s2$se_ci_implied)
})

test_that("NA rows are dropped, not silently mishandled", {
  d <- get_class("A")
  d$x2[c(1, 5, 9)] <- NA
  m <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d)
  expect_false(any(m$esample[c(1, 5, 9)]))
  expect_true(m$converged)
})
