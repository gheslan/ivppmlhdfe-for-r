# Multi-part formula parsing.

test_that("FE and IV parts may appear in either order", {
  d <- get_class("A")
  m1 <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d)
  m2 <- ivppmlhdfe(y ~ x2 | x1 ~ z | id + year, data = d)
  expect_equal(coef(m1), coef(m2), tolerance = 1e-10)
})

test_that("no exogenous regressors (y ~ 1 | fe | endog ~ z)", {
  d <- get_class("A")
  m <- ivppmlhdfe(y ~ 1 | id + year | x1 ~ z, data = d)
  expect_setequal(names(coef(m)), c("x1", "_cons"))
})

test_that("no fixed effects (Stata noabsorb) works", {
  d <- get_class("A")
  m <- ivppmlhdfe(y ~ x2 | x1 ~ z, data = d)
  expect_true(m$converged)
  expect_equal(m$df_a, 0L)
})

test_that("interaction FE via ^ equals manual interaction id", {
  d <- get_class("B")
  m1 <- ivppmlhdfe(y ~ x2 | exp^year + imp^year | x1 ~ z, data = d)
  d$ey <- interaction(d$exp, d$year, drop = TRUE)
  d$iy <- interaction(d$imp, d$year, drop = TRUE)
  m2 <- ivppmlhdfe(y ~ x2 | ey + iy | x1 ~ z, data = d)
  expect_equal(unname(coef(m1)["x1"]), unname(coef(m2)["x1"]), tolerance = 1e-8)
})

test_that("# is accepted as interaction operator (Stata syntax)", {
  d <- get_class("B")
  m1 <- ivppmlhdfe(y ~ x2 | exp^year + imp^year | x1 ~ z, data = d)
  m2 <- ivppmlhdfe(y ~ x2 | `exp#year` + `imp#year` | x1 ~ z, data = d)
  expect_equal(unname(coef(m1)["x1"]), unname(coef(m2)["x1"]), tolerance = 1e-10)
})

test_that("too many | parts is rejected", {
  d <- get_class("A")
  expect_error(ivppmlhdfe(y ~ x2 | id | year | x1 ~ z, data = d),
               "too many")
})

test_that("subset (Stata if/in) restricts the sample", {
  d <- get_class("A")
  m <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d,
                  subset = d$year <= 5)
  # esample is defined over the subsetted data frame
  expect_lte(length(m$esample), sum(d$year <= 5))
  expect_lte(m$N, sum(d$year <= 5))
})
