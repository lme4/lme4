library("testthat")
library("lme4")

L <- load(system.file("testdata", "lme-tst-fits.rda",
                      package="lme4", mustWork=TRUE))

context("Show basic allFit results")
test_that("allFit print/summary is fine", {
    gm_all <- allFit(fit_cbpp_1, verbose=FALSE)
    expect_is(gm_all, "allFit")
    expect_is(summary(gm_all), "summary.allFit")
})
