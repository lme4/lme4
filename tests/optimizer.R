library(lme4)
(testLevel <- if (nzchar(s <- Sys.getenv("LME4_TEST_LEVEL"))) as.numeric(s) else 1)
source(system.file("test-tools-1.R", package = "Matrix"), keep.source = FALSE)

## should be able to run any example with any bounds-constrained optimizer ...
##  Nelder_Mead, bobyqa built in; optimx/nlminb, optimx/L-BFGS-B
##  optimx/Rcgmin will require a bit more wrapping/interface work (requires gradient)

fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)  ## Nelder_Mead
fm1B <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy,
             control=lmerControl(optimizer="bobyqa"))
stopifnot(all.equal(fixef(fm1),fixef(fm1B)))
require(optimx)
lmerCtrl.optx <- function(method, ...)
    lmerControl(optimizer="optimx", ..., optCtrl=list(method=method))
glmerCtrl.optx <- function(method, ...)
    glmerControl(optimizer="optimx", ..., optCtrl=list(method=method))

## FAILS on Windows (on r-forge only, not win-builder)... 'function is infeasible at initial parameters'
## (can we test whether we are on r-forge??)
if (.Platform$OS.type != "windows") {
    fm1C <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy,
                 control=lmerCtrl.optx(method="nlminb"))
    fm1D <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy,
                 control=lmerCtrl.optx(method="L-BFGS-B"))
    stopifnot(is.all.equal4(fixef(fm1),fixef(fm1B),fixef(fm1C),fixef(fm1D)))

    if (testLevel > 2) {
        fm1E <- update(fm1,control=lmerCtrl.optx(method=c("nlminb","L-BFGS-B")))
        ## hack equivalence of call and optinfo
        fm1E@call <- fm1C@call
        fm1E@optinfo <- fm1C@optinfo
        ## FIXME: this *should* be identical, but we have small numeric differences
        all.equal(fm1C,fm1E,tol=1e-5)
    }
}

gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
             data = cbpp, family = binomial,
             control=glmerControl(tolPwrss=1e-13))
gm1B <- update(gm1, control=glmerControl  (tolPwrss=1e-13, optimizer="bobyqa"))
gm1C <- update(gm1, control=glmerCtrl.optx(tolPwrss=1e-13, method="nlminb"))
gm1D <- update(gm1, control=glmerCtrl.optx(tolPwrss=1e-13, method="L-BFGS-B"))
stopifnot(is.all.equal4(fixef(gm1),fixef(gm1B),fixef(gm1C),fixef(gm1D),tol=1e-5))

gm1E <- update(gm1, control=glmerCtrl.optx(tolPwrss=1e-13,
                    method=c("nlminb","L-BFGS-B")))


if (testLevel > 2) {
    ## hack equivalence of call and optinfo
    gm1E@call <- gm1C@call
    gm1E@optinfo <- gm1C@optinfo
    all.equal(gm1E,gm1C,tol=1e-5)  ## FIXME: this *should* be identical, but we have small numeric differences
}
