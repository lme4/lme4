library("testthat")
library("lme4")

context("storing warnings, convergence status, etc.")

test_that("storewarning", {
    expect_warning(glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
          cbpp, binomial, control=list(maxfun=3),
          optimizer="Nelder_Mead"),"failure to converge in 3")
})
