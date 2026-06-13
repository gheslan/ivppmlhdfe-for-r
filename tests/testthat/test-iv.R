# Instrumental-variable behaviour.

test_that("IV removes endogeneity bias (Class A, true beta1 = 0.5)", {
  d <- get_class("A")
  iv <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d)
  ppml <- ivppmlhdfe(y ~ x1 + x2 | id + year, data = d)   # no IV -> biased
  # IV estimate close to truth; PPML badly biased upward
  expect_lt(abs(coef(iv)["x1"] - 0.5), 0.05)
  expect_gt(coef(ppml)["x1"], 0.6)
})

test_that("just-identified IV runs on all three classes", {
  for (cl in c("A", "B", "C")) {
    d <- get_class(cl)
    fml <- switch(cl,
      A = y ~ x2 | id + year | x1 ~ z,
      B = y ~ x2 | exp^year + imp^year | x1 ~ z,
      C = y ~ x2 | exp^year + imp^year + pair | x1 ~ z)
    m <- ivppmlhdfe(fml, data = d)
    expect_true(m$converged)
    expect_true("x1" %in% m$endog_names)
    expect_true("z" %in% m$inst_names)
  }
})

test_that("under-identification is rejected", {
  d <- get_class("A")
  d$z2 <- d$z + d$x2
  # two endogenous, one instrument -> not identified
  expect_error(
    ivppml_fit(y = d$y, X_endog = cbind(d$x1, d$x2), Z_excl = matrix(d$z),
               fe = list(d$id, d$year)),
    "not identified")
})

test_that("instrument equal to a regressor is rejected", {
  d <- get_class("A")
  expect_error(
    ivppml_fit(y = d$y, X_exog = matrix(d$x2, ncol = 1),
               X_endog = matrix(d$x1, ncol = 1),
               Z_excl = matrix(d$x1, ncol = 1),  # z == endog regressor
               fe = list(d$id, d$year)),
    "instrument cannot equal a regressor")
})

test_that("matrix and formula interfaces agree", {
  d <- get_class("A")
  mf <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d)
  mm <- ivppml_fit(y = d$y, X_exog = matrix(d$x2, ncol = 1,
                                            dimnames = list(NULL, "x2")),
                   X_endog = matrix(d$x1, ncol = 1,
                                    dimnames = list(NULL, "x1")),
                   Z_excl = matrix(d$z, ncol = 1, dimnames = list(NULL, "z")),
                   fe = list(id = d$id, year = d$year))
  expect_equal(unname(coef(mf)[c("x2", "x1")]),
               unname(coef(mm)[c("x2", "x1")]), tolerance = 1e-10)
})
