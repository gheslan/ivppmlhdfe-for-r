# Split-panel jackknife bias correction.

test_that("Class A SPJ formula uses 3*b_full - mean(T) - mean(N)", {
  d <- get_class("A")
  s <- ivppml_spj(y ~ x2 | id + year | x1 ~ z, data = d,
                  class = "A", id = "id", time = "year", B = 0, seed = 1)
  sub <- s$subpanel
  expected <- 3 * s$b_full - (sub$T1 + sub$T2) / 2 - (sub$N1 + sub$N2) / 2
  expect_equal(unname(s$b_spj), unname(expected), tolerance = 1e-12)
  expect_s3_class(s, "ivppml_spj")
})

test_that("Class B SPJ formula uses 2*b_full - mean(4 country panels)", {
  d <- get_class("B")
  s <- ivppml_spj(y ~ x2 | exp^year + imp^year | x1 ~ z, data = d,
                  class = "B", exp = "exp", imp = "imp", B = 0, seed = 2)
  expected <- 2 * s$b_full - mean(s$subpanel$country)
  expect_equal(unname(s$b_spj), unname(expected), tolerance = 1e-12)
})

test_that("Class C SPJ formula uses 4*b - 2*country - 2*time + 8cell", {
  d <- get_class("C")
  s <- ivppml_spj(y ~ x2 | exp^year + imp^year + pair | x1 ~ z, data = d,
                  class = "C", exp = "exp", imp = "imp", time = "year",
                  B = 0, seed = 3)
  sp <- s$subpanel
  expected <- 4 * s$b_full - 2 * mean(sp$country) - 2 * mean(sp$time) +
    mean(sp$cell8)
  expect_equal(unname(s$b_spj), unname(expected), tolerance = 1e-12)
})

test_that("interacted single-family SPJ uses 2*b_full - mean(N)", {
  d <- get_class("A")
  # single interacted FE family (year) that pools across the split unit (id);
  # the cross-sectional split is on id, as in the Chen et al. (2025) replication
  s <- ivppml_spj(y ~ x2 | year | x1 ~ z, data = d,
                  class = "interacted", id = "id", B = 0, seed = 4)
  expected <- 2 * s$b_full - (s$subpanel$N1 + s$subpanel$N2) / 2
  expect_equal(unname(s$b_spj), unname(expected), tolerance = 1e-12)
})

test_that("SPJ requires the relevant id/time/exp/imp arguments", {
  d <- get_class("A")
  expect_error(ivppml_spj(y ~ x2 | id + year | x1 ~ z, data = d, class = "A",
                          time = "year"), "requires `id`")
})

test_that("bootstrap produces a CI-implied SE and percentile CI", {
  d <- get_class("A")
  s <- ivppml_spj(y ~ x2 | id + year | x1 ~ z, data = d,
                  class = "A", id = "id", time = "year", B = 8, seed = 5)
  expect_true(is.finite(s$se_ci_implied) || s$n_failed == s$B)
  expect_length(s$ci, 2L)
})
