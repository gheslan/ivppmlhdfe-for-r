# The sandwich (bread, meat) must equal fepois exactly; only the small-sample
# multiplier differs (ivppmlhdfe uses N/(N-1) robust, G/(G-1) cluster).

test_that("robust sandwich equals fepois up to N/(N-1) multiplier", {
  skip_no_fixest()
  d <- get_class("A")
  m <- ivppmlhdfe(y ~ x1 + x2 | id + year, data = d)
  fe <- fixest::fepois(y ~ x1 + x2 | id + year, data = d, notes = FALSE)
  N <- m$N
  se_ours <- sqrt(diag(vcov(m)))[c("x1", "x2")]
  # fepois pure sandwich (no small-sample adjustment) -> multiplier 1
  se_fe <- sqrt(diag(vcov(fe, vcov = "hetero",
                                  ssc = fixest::ssc(adj = FALSE,
                                                    cluster.adj = FALSE))))[c("x1", "x2")]
  ratio <- unname(se_ours / se_fe)
  expect_equal(ratio, rep(sqrt(N / (N - 1)), 2), tolerance = 1e-6)
})

test_that("cluster-robust VCE is bit-identical to fepois (G/(G-1))", {
  skip_no_fixest()
  d <- get_class("B")
  m <- ivppmlhdfe(y ~ x1 + x2 | exp^year + imp^year, data = d, cluster = ~pair)
  fe <- fixest::fepois(y ~ x1 + x2 | exp^year + imp^year, data = d, notes = FALSE)
  se_o <- sqrt(diag(vcov(m)))[c("x1", "x2")]
  se_f <- sqrt(diag(vcov(fe, cluster = ~pair,
                                 ssc = fixest::ssc(adj = FALSE,
                                                   cluster.adj = TRUE))))[c("x1", "x2")]
  expect_equal(unname(se_o), unname(se_f), tolerance = 1e-6)
})

test_that("two-way clustering runs and is symmetric PSD", {
  d <- get_class("B")
  m <- ivppmlhdfe(y ~ x2 | exp^year + imp^year | x1 ~ z, data = d,
                  cluster = ~exp + imp)
  V <- vcov(m)[c("x1", "x2"), c("x1", "x2"), drop = FALSE]
  expect_equal(V, t(V), tolerance = 1e-12)
  expect_true(min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) >= -1e-10)
})

test_that("single cluster is demoted to robust with a warning", {
  d <- get_class("A")
  d$onecl <- 1L
  expect_warning(
    m <- ivppmlhdfe(y ~ x1 + x2 | id + year, data = d, cluster = ~onecl),
    "1 cluster"
  )
  expect_equal(m$vcetype, "robust")
})
