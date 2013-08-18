library("testthat")
library("lme4")

context("testing bounds")

test_that("lmerBounds", {
    lmod <- lFormula(Reaction ~ Days + (Days|Subject), sleepstudy)
    devfun <- do.call(mkLmerDevfun, lmod)
    opt1 <- optimizeLmer(devfun)
    environment(devfun)$pp$theta <- opt1$par
    m1 <- mkMerMod(environment(devfun), opt1, lmod$reTrms, fr = lmod$fr)
    ## 
    devfun2 <- do.call(mkLmerDevfun, lmod)
    environment(devfun2)$upper <- c(1.1,0.005,1.1)
    opt2 <- optimizeLmer(devfun2)
    environment(devfun2)$pp$theta <- opt2$par
    m2 <- mkMerMod(environment(devfun2), opt2, lmod$reTrms, fr = lmod$fr)
})



