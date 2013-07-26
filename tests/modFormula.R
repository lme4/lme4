library(lme4)
library(testthat)
lmod <- lFormula(Reaction ~ Days + (Days|Subject), sleepstudy)
devfun <- do.call(mkLmerDevfun, lmod)
opt <- optimizeLmer(devfun)
fm1 <- mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr)
fm2 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)

## basic equivalence
fm1C <- fm1
fm1C@call <- fm2@call
expect_equal(fm2,fm1C)
expect_equal(range(residuals(fm1)),c(-101.1789,132.5466),tol=1e-6)
expect_is(model.frame(fm1),"data.frame")
## formulae
expect_equal(formula(model.frame(fm1)),Reaction ~ Days + Subject) ## fixed only
expect_equal(formula(fm1),Reaction~Days+(Days|Subject))
## predictions
expect_equal(predict(fm1,newdata=sleepstudy[1:10,],REform=NULL),
             predict(fm2,newdata=sleepstudy[1:10,],REform=NULL))
expect_equal(unname(predict(fm1,newdata=sleepstudy)),
             predict(fm1))

lmodOff <- lFormula(Reaction ~ Days + (Days|Subject) + offset(0.5*Days),
                 sleepstudy)
devfunOff <- do.call(mkLmerDevfun, lmodOff)
opt <- optimizeLmer(devfunOff)
fm1Off <- mkMerMod(environment(devfunOff), opt, lmodOff$reTrms, fr = lmodOff$fr)
fm2Off <- lmer(Reaction ~ Days + (Days|Subject) + offset(0.5*Days), sleepstudy)
expect_equal(predict(fm1Off,newdata=sleepstudy[1:10,],REform=NULL),
             predict(fm2Off,newdata=sleepstudy[1:10,],REform=NULL))

## FIXME: need more torture tests with offset specified, in different environments ...

## FIXME: drop1(.) doesn't work with modular objects ... hard to see how it
##  could, though ...
## drop1(fm1Off)
drop1(fm2Off)
