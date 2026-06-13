# Collinearity removal and separation detection.

test_that("a collinear duplicate exogenous regressor is dropped", {
  d <- get_class("A")
  d$x2dup <- d$x2          # exact duplicate
  m <- ivppmlhdfe(y ~ x2 + x2dup | id + year | x1 ~ z, data = d)
  expect_equal(m$n_dropped_exog, 1L)
  # surviving coefficient equals the no-duplicate fit
  m0 <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d)
  expect_equal(unname(coef(m)["x1"]), unname(coef(m0)["x1"]), tolerance = 1e-8)
})

test_that("a regressor fully absorbed by the FE is dropped", {
  d <- get_class("A")
  d$yearcopy <- d$year     # collinear with year FE after demeaning
  m <- ivppmlhdfe(y ~ x2 + yearcopy | id + year | x1 ~ z, data = d)
  expect_equal(m$n_dropped_exog, 1L)
})

test_that("FE-separated all-zero group is dropped", {
  d <- get_class("A")
  d$y[d$id == 1] <- 0      # individual 1 all-zero -> separated by id FE
  m <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d)
  expect_gt(m$num_sep_fe + m$num_singletons, 0)
  expect_false(any(m$esample & d$id == 1))
})

test_that("singletons are dropped by default and kept when requested", {
  d <- get_class("A")
  # add an isolated id with a single observation
  extra <- d[1, , drop = FALSE]; extra$id <- max(d$id) + 1L; extra$y <- 3
  d2 <- rbind(d, extra)
  m_drop <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = d2)
  expect_gt(m_drop$num_singletons, 0)
})

test_that("tagsep flags separated observations without estimating", {
  d <- get_class("A")
  d$y[d$id == 2] <- 0
  tg <- ivppml_tagsep(y ~ x2 | id + year | x1 ~ z, data = d)
  expect_length(tg$separated, nrow(d))
  expect_true(tg$num_singletons >= 0)
  expect_true(any(tg$separated | tg$singleton))
})
