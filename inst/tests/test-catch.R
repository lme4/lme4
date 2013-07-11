library("testthat")
library("lme4")

context("storing warnings, convergence status, etc.")

test_that("storewarning", {
    expect_warning(gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                         data=cbpp, family=binomial,
                         control=glmerControl(optimizer="Nelder_Mead",maxfun=3)),
                         "failure to converge in 3")
    expect_equal(gm1@optinfo$warnings[[1]],"failure to converge in 3 evaluations")
    ## FIXME: why is conv==0 here?
})
