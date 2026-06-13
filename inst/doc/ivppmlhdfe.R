## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## -----------------------------------------------------------------------------
library(ivppmlhdfe)
data(ivppmlhdfe_ClassA)

## -----------------------------------------------------------------------------
ppml <- ivppmlhdfe(y ~ x1 + x2 | id + year, data = ivppmlhdfe_ClassA)
coef(ppml)["x1"]   # biased upward, away from 0.5

## -----------------------------------------------------------------------------
iv <- ivppmlhdfe(y ~ x2 | id + year | x1 ~ z, data = ivppmlhdfe_ClassA)
coef(iv)["x1"]     # close to the true 0.5
summary(iv)

## -----------------------------------------------------------------------------
data(ivppmlhdfe_ClassB)
b <- ivppmlhdfe(y ~ x2 | exp^year + imp^year | x1 ~ z,
                data = ivppmlhdfe_ClassB, cluster = ~pair)
coef(b)["x1"]

## -----------------------------------------------------------------------------
s <- ivppml_spj(y ~ x2 | id + year | x1 ~ z, data = ivppmlhdfe_ClassA,
                class = "A", id = "id", time = "year", B = 0)
s$b_full   # uncorrected
s$b_spj    # SPJ-corrected

## -----------------------------------------------------------------------------
if (requireNamespace("fixest", quietly = TRUE)) {
  fe <- fixest::fepois(y ~ x1 + x2 | id + year, data = ivppmlhdfe_ClassA,
                       notes = FALSE)
  rbind(ivppmlhdfe = coef(ppml)[c("x1", "x2")],
        fepois     = coef(fe)[c("x1", "x2")])
}

