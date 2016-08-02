library("testthat")
library("lme4")

context("storing warnings, convergence status, etc.")

test_that("storewarning", {
    gCtrl <- glmerControl(optimizer = "Nelder_Mead",
                          optCtrl = list(maxfun=3))
    expect_warning(gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                                data=cbpp, family=binomial, control=gCtrl),
                   "failure to converge in 3")
    expect_equal(gm1@optinfo$warnings[[1]],"failure to converge in 3 evaluations")
    ## conv==0 here because Nelder_Mead doesn't set conv ...
})
