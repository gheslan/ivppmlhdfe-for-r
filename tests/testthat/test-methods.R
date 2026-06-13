# S3 methods: coef, vcov, predict, summary, confint, residuals, logLik, nobs.

setup_fit <- function() {
  d <- get_class("A")
  ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d)
}

test_that("accessors return correct shapes", {
  m <- setup_fit()
  expect_named(coef(m), c("x2", "x1", "_cons"))
  expect_equal(dim(vcov(m)), c(3L, 3L))
  expect_equal(nobs(m), m$N)
  expect_s3_class(logLik(m), "logLik")
})

test_that("predict types are consistent", {
  m <- setup_fit()
  mu  <- predict(m, "mu")
  xbd <- predict(m, "xbd")
  d   <- predict(m, "d")
  xb  <- predict(m, "xb")
  expect_equal(mu, exp(xbd), tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(xbd, xb + d, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(predict(m, "response"), m$y - mu, tolerance = 1e-12, ignore_attr = TRUE)
  # mu equals stored fitted values
  expect_equal(mu, fitted(m), tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("predict stdp is nonnegative and finite", {
  m <- setup_fit()
  s <- predict(m, "stdp")
  expect_length(s, m$N)
  expect_true(all(is.finite(s)))
  expect_true(all(s >= 0))
})

test_that("predict deviance contributions reproduce total deviance", {
  m <- setup_fit()
  dc <- predict(m, "deviance")
  expect_equal(sum(dc), m$deviance, tolerance = 1e-6)
})

test_that("residual types are internally consistent", {
  m <- setup_fit()
  mu <- fitted(m); y <- m$y
  expect_equal(residuals(m, "response"), y - mu, tolerance = 1e-12)
  expect_equal(residuals(m, "pearson"), (y - mu) / sqrt(mu), tolerance = 1e-12)
  expect_equal(residuals(m, "working"), (y - mu) / mu, tolerance = 1e-12)
})

test_that("confint width matches z * se", {
  m <- setup_fit()
  ci <- confint(m)
  se <- sqrt(diag(vcov(m)))
  half <- (ci[, 2] - ci[, 1]) / 2
  expect_equal(unname(half), unname(stats::qnorm(0.975) * se), tolerance = 1e-10)
})

test_that("summary builds a coefficient table; eform exponentiates", {
  m <- setup_fit()
  s <- summary(m)
  expect_s3_class(s, "summary.ivppmlhdfe")
  expect_true(all(c("Coef.", "Std. Err.", "z", "P>|z|") %in% colnames(s$coeftable)))
  se <- summary(m, eform = TRUE)
  expect_equal(unname(se$coeftable[, 1]), unname(exp(coef(m))), tolerance = 1e-10)
})

test_that("print methods run without error", {
  m <- setup_fit()
  expect_output(print(m))
  expect_output(print(summary(m)))
})
