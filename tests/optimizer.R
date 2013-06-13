library(lme4)
source(system.file("test-tools-1.R", package = "Matrix"), keep.source = FALSE)

## should be able to run any example with any bounds-constrained optimizer ...
##  Nelder_Mead, bobyqa built in; optimx/nlminb, optimx/L-BFGS-B
##  optimx/Rcgmin will require a bit more wrapping/interface work (requires gradient)

fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)  ## Nelder_Mead
fm1B <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy,
             control=lmerControl(optimizer="bobyqa"))
stopifnot(all.equal(fixef(fm1),fixef(fm1B)))
require(optimx)
## FAILS on Windows (on r-forge only, not win-builder)... 'function is infeasible at initial parameters'
## (can we test whether we are on r-forge??)
if (.Platform$OS.type != "windows") {
    fm1C <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy,
                 control=lmerControl(optimizer="optimx",method="nlminb"))
    fm1D <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy,
                 control=lmerControl(optimizer="optimx", method="L-BFGS-B"))
    stopifnot(is.all.equal4(fixef(fm1),fixef(fm1B),fixef(fm1C),fixef(fm1D)))
}

gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
             data = cbpp, family = binomial,
             control=glmerControl(tolPwrss=1e-13))
gm1B <- update(gm1,control=glmerControl(tolPwrss=1e-13,optimizer="bobyqa"))
gm1C <- update(gm1,control=glmerControl(optimizer="optimx",method="nlminb",tolPwrss=1e-13))
gm1D <- update(gm1,control=glmerControl(optimizer="optimx",method="L-BFGS-B",tolPwrss=1e-13))
stopifnot(is.all.equal4(fixef(gm1),fixef(gm1B),fixef(gm1C),fixef(gm1D),tol=1e-5))


