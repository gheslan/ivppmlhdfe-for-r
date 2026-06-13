# Plain PPML (no IV) must reproduce ppmlhdfe / fixest::fepois exactly.
# This is the strongest available external numerical benchmark.

test_that("plain PPML coefficients match fepois to machine precision (Class A)", {
  skip_no_fixest()
  d <- get_class("A")
  m <- ivppmlhdfe(y ~ x1 + x2 | id + year, data = d)
  fe <- fixest::fepois(y ~ x1 + x2 | id + year, data = d, notes = FALSE)
  cb <- coef(m)[c("x1", "x2")]
  fb <- coef(fe)[c("x1", "x2")]
  expect_equal(unname(cb), unname(fb), tolerance = 1e-10)
})

test_that("plain PPML matches fepois for Class B (two-way gravity)", {
  skip_no_fixest()
  d <- get_class("B")
  m <- ivppmlhdfe(y ~ x1 + x2 | exp^year + imp^year, data = d)
  fe <- fixest::fepois(y ~ x1 + x2 | exp^year + imp^year, data = d, notes = FALSE)
  expect_equal(unname(coef(m)[c("x1", "x2")]),
               unname(coef(fe)[c("x1", "x2")]), tolerance = 1e-9)
})

test_that("plain PPML matches fepois for Class C (three-way gravity)", {
  skip_no_fixest()
  d <- get_class("C")
  m <- ivppmlhdfe(y ~ x1 + x2 | exp^year + imp^year + pair, data = d)
  fe <- fixest::fepois(y ~ x1 + x2 | exp^year + imp^year + pair, data = d,
                       notes = FALSE)
  expect_equal(unname(coef(m)[c("x1", "x2")]),
               unname(coef(fe)[c("x1", "x2")]), tolerance = 1e-7)
})

test_that("log-likelihood matches fepois", {
  skip_no_fixest()
  d <- get_class("A")
  m <- ivppmlhdfe(y ~ x1 + x2 | id + year, data = d)
  fe <- fixest::fepois(y ~ x1 + x2 | id + year, data = d, notes = FALSE)
  expect_equal(m$ll, as.numeric(stats::logLik(fe)), tolerance = 1e-8)
})

test_that("constant is recovered (ivreg2 partial convention, zero variance)", {
  d <- get_class("A")
  m <- ivppmlhdfe(y ~ x1 + x2 | id + year, data = d)
  expect_true("_cons" %in% names(coef(m)))
  expect_equal(unname(diag(vcov(m))["_cons"]), 0)
})
