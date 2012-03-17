library(lme4)
## should be able to run any example with any bounds-constrained optimizer ...
##  Nelder_Mead, bobyqa built in; optimx/nlminb, optimx/L-BFGS-B
##  optimx/Rcgmin will require a bit more wrapping/interface work (requires gradient)

fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)  ## Nelder_Mead
fm1B <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy,
            optimizer="bobyqa")
stopifnot(all.equal(fixef(fm1),fixef(fm1B)))
require(optimx)
## FAILS on Windows (on r-forge only, not win-builder)... 'function is infeasible at initial parameters'
## (can we test whether we are on r-forge??)
if (.Platform$OS.type != "windows") {
    fm1C <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy,
                 optimizer="optimx", control=list(method="nlminb"))
    fm1D <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy,
                 optimizer="optimx", control=list(method="L-BFGS-B"))
    stopifnot(all.equal(fixef(fm1),fixef(fm1B),fixef(fm1C),fixef(fm1D)))
}

gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                   data = cbpp, family = binomial, tolPwrss=1e-13)
gm1B <- update(gm1,optimizer="bobyqa")

if (FALSE) {
  gm1C <- update(gm1,optimizer="optimx",control=list(method="nlminb"))
  gm1D <- update(gm1,optimizer="optimx",control=list(method="L-BFGS-B"))

  ##  these fail because the initial optimization step finds the (true)
  ## maximum of the nAGQ-0 deviance function at zero; then nlminb and L-BFGS-B
  ## can't get off the boundary ... don't know how well they would do even if
  ## they could

  gm1fun <- update(gm1,devFunOnly=TRUE)
  gm1fun0 <- update(gm1,devFunOnly=TRUE,nAGQ=0)
  svec <- seq(0.01,3,by=0.01)
  dvec <- sapply(svec,gm1fun0)
  plot(svec,dvec,type="l")
  abline(v=1)
  n1 <- nlminb(start=1,objective=gm1fun0,control=list(),lower=0)  ## false convergence=code 1
  b1 <- bobyqa(par=1,fn=gm1fun0,control=list(),lower=0)
}

