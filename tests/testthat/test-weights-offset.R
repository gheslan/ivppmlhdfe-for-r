# Weights, offset and exposure semantics.

test_that("frequency weights set effective N = sum(weights)", {
  d <- get_class("A")
  w <- rep(2, nrow(d))
  m <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d,
                  weights = w, weight_type = "fweight")
  expect_equal(m$N, sum(w))
})

test_that("unit pweights reproduce the unweighted fit", {
  d <- get_class("A")
  m0 <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d)
  m1 <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d,
                   weights = rep(1, nrow(d)))
  expect_equal(coef(m0), coef(m1), tolerance = 1e-10)
})

test_that("pweights change estimates and agree with fepois point estimates", {
  skip_no_fixest()
  d <- get_class("A")
  set.seed(7); w <- runif(nrow(d), 0.5, 2)
  m <- ivppmlhdfe(y ~ x1 + x2 | id + year, data = d, weights = w)
  fe <- fixest::fepois(y ~ x1 + x2 | id + year, data = d, weights = w,
                       notes = FALSE)
  expect_equal(unname(coef(m)[c("x1", "x2")]),
               unname(coef(fe)[c("x1", "x2")]), tolerance = 1e-8)
})

test_that("offset with coefficient 1 matches fepois offset", {
  skip_no_fixest()
  d <- get_class("A")
  set.seed(11); off <- log(runif(nrow(d), 1, 3))
  m <- ivppmlhdfe(y ~ x1 + x2 | id + year, data = d, offset = off)
  fe <- fixest::fepois(y ~ x1 + x2 | id + year, data = d,
                       offset = off, notes = FALSE)
  expect_equal(unname(coef(m)[c("x1", "x2")]),
               unname(coef(fe)[c("x1", "x2")]), tolerance = 1e-7)
})

test_that("exposure(v) equals offset(log(v))", {
  d <- get_class("A")
  set.seed(13); v <- runif(nrow(d), 1, 4)
  m_e <- ivppmlhdfe(y ~ x1 + x2 | id + year, data = d, exposure = v)
  m_o <- ivppmlhdfe(y ~ x1 + x2 | id + year, data = d, offset = log(v))
  expect_equal(coef(m_e), coef(m_o), tolerance = 1e-10)
})
