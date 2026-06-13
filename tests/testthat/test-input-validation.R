# Input validation and informative errors.

test_that("negative outcomes are rejected", {
  d <- get_class("A"); d$y[1] <- -1
  expect_error(ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d),
               "nonnegative")
})

test_that("all-zero outcome is rejected", {
  d <- get_class("A"); d$y <- 0
  expect_error(ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d),
               "identically zero")
})

test_that("non-finite regressors are rejected", {
  d <- get_class("A"); d$x2[5] <- Inf
  expect_error(ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d))
})

test_that("constant instrument is rejected", {
  d <- get_class("A")
  expect_error(
    ivppml_fit(y = d$y, X_endog = matrix(d$x1, ncol = 1),
               Z_excl = matrix(1, nrow(d), 1), fe = list(d$id, d$year)),
    "zero variance")
})

test_that("negative weights are rejected", {
  d <- get_class("A"); w <- rep(1, nrow(d)); w[1] <- -2
  expect_error(ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d, weights = w),
               "negative weights")
})

test_that("exposure must be positive; offset and exposure are mutually exclusive", {
  d <- get_class("A")
  expect_error(ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d,
                          exposure = rep(-1, nrow(d))), "greater than zero")
  expect_error(ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d,
                          exposure = rep(1, nrow(d)), offset = rep(0, nrow(d))),
               "both")
})

test_that("unknown separation technique is rejected", {
  d <- get_class("A")
  expect_error(ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d,
                          separation = "banana"), "unknown technique")
})

test_that("non-data.frame data is rejected", {
  expect_error(ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = 1:10),
               "data.frame")
})
