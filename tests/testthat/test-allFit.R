library("testthat")
library("lme4")

L <- load(system.file("testdata", "lme-tst-fits.rda",
                      package="lme4", mustWork=TRUE))

gm_all <- allFit(fit_cbpp_1, verbose=FALSE)

context("Show basic allFit results")
test_that("allFit print/summary is fine", {
    expect_is(gm_all, "allFit")
    expect_is(summary(gm_all), "summary.allFit")
})

test_that("nloptwrap switches optimizer correctly", {
    expect_equal(attr(gm_all[["nloptwrap.NLOPT_LN_BOBYQA"]],"optCtrl"),
                 list(algorithm = "NLOPT_LN_BOBYQA"))
    expect_equal(attr(gm_all[["nloptwrap.NLOPT_LN_NELDERMEAD"]],"optCtrl"),
                 list(algorithm = "NLOPT_LN_NELDERMEAD"))

})

test_that("lmerControl() arg works too", {
    fm0 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
    fm  <- update(fm0, control = lmerControl(optCtrl = list(xtol_rel = 1e-8,
                                                            ftol_rel = 1e-12)))
    afm0 <- allFit(fm0,verbose=FALSE)
    afm  <- allFit(fm,verbose=FALSE) # used to fail
    drop_ <- function(x) {
        x[setdiff(names(x), c("times","feval"))]
    }
    ## print(all.equal(drop_(summary(afm0)), drop_(summary(afm)), tolerance = 0))
    expect_equal(drop_(summary(afm0)),
                 drop_(summary(afm)), tolerance = 1e-3)
})

test_that("glmerControl() arg + optimizer", {
    ## GH #523?
    fit_cbpp_1u <- update(fit_cbpp_1,
                          control=glmerControl(optimizer="nloptwrap",
                         optCtrl=list(xtol_abs=1e-10, ftol_abs=1e-10)))
    af2 <- allFit(fit_cbpp_1u, verbose=FALSE)
    expect_equal(class(af2),"allFit")
})
