library("testthat")
library("lme4")
context("specifying starting values")

##' Update 'mod', copying .@call and attr(.@frame, "start")  from 'from'
copysome <- function(mod, from) {
    stopifnot(all.equal(class(mod), class(from)), isS4(mod))
    mod@call <- from@call
    attr(mod@frame, "start") <- attr(from@frame, "start")
    mod
}

## is "Nelder_Mead" default optimizer?
isNM <- formals(lmerControl)$optimizer == "Nelder_Mead"

stMsg <- "'start' must be .* a numeric vector .* list"
test_that("lmer", {
    frm <- Reaction ~ Days + (Days|Subject)
    ctrl <- lmerControl(optCtrl = list(maxfun= if(isNM) 50 else 100))
    x <- suppressWarnings(lmer(frm, data=sleepstudy, control=ctrl, REML=FALSE))
    x2 <- suppressWarnings(update(x,start=c(1,0,1)))
    x3 <- suppressWarnings(update(x,start=list(theta=c(1,0,1))))
    ff <- update(x,devFunOnly=TRUE)
    x2@call <- x3@call <- x@call  ## hack call component
    expect_equal(x,x2)
    expect_equal(x,x3)
    expect_error(update(x, start = "a"), stMsg)
    ## misspelled
    expect_error(update(x,start=list(Theta=c(1,0,1))),"incorrect components")
    th0 <- getME(x,"theta")
    y <- suppressWarnings(update(x,start=th0))
    if(isNM) {
        expect_equal(AIC(x), 1768.025, tolerance=1e-6)
        expect_equal(AIC(y), 1763.949, tolerance=1e-6)
    } else { ## only "bobyqa" tested:
        expect_equal(AIC(x), 1763.939344, tolerance=1e-6)
        expect_equal(AIC(x), AIC(y))
    }
    if(isNM)
        expect_equal(suppressWarnings(optimizeLmer(ff,control=list(maxfun=50),
                                                   start=c(1,0,1))$fval),
                 unname(deviance(x)))
    expect_equal(suppressWarnings(optimizeLmer(ff,control=list(maxfun=50),
                                               start=th0)$fval),
                 unname(deviance(y)))
})
test_that("glmer", {
    ctrl <- glmerControl(optCtrl=list(maxfun=50)) # -> non-convergence warnings
    x <- suppressWarnings(glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                   data = cbpp, family = binomial, control=ctrl))
    ## theta only
    x2 <- suppressWarnings(update(x, start= 1))
    x3 <- suppressWarnings(update(x, start= list(theta = 1)))
    ff <- update(x,devFunOnly=TRUE)
    x2@call <- x3@call <- x@call  ## hack call component
    expect_equal(x,x2)
    expect_equal(x,x3)
    expect_error(update(x, start="a"), stMsg)
    expect_error(update(x, start=list(Theta=1)), "bad name\\(s\\)")
    th0 <- getME(x,"theta")
    y <- suppressWarnings(update(x, start=th0)) # expect_equal() fails: optinfo -> derivs -> Hessian

    ## theta and beta
    x0 <- update(x,nAGQ=0)
    x4 <- suppressWarnings(update(x, start = list(theta=1, fixef=fixef(x0))))
    expect_equal(x, copysome(x4, from=x))
    x5 <- suppressWarnings(update(x, start = list(theta=1, fixef=rep(0,4))))
    expect_equal(AIC(x5), 221.5823, tolerance=1e-6)
    x6 <- expect_error(update(x, start = list(theta=1, fixef=rep(0,5))),
                       "incorrect number of fixef components")
    ## beta only
    x7 <- suppressWarnings(update(x,start=list(fixef=fixef(x0))))
    expect_equal(x, copysome(x7, from=x))
    x8 <- suppressWarnings(update(x,start=list(fixef=rep(0,4))))
    expect_equal(x5, copysome(x8, from=x5))
})
