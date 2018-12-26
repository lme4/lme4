library(lme4)
library(testthat)

.get.checkingOpts <- lme4:::.get.checkingOpts
stopifnot(identical(
    .get.checkingOpts(
    c("CheckMe", "check.foo", "check.conv.1", "check.rankZ", "check.rankX"))
    , c("check.foo", "check.rankZ")))

lmod <- lFormula(Reaction ~ Days + (Days|Subject), sleepstudy)
devfun <- do.call(mkLmerDevfun, lmod)
opt <- optimizeLmer(devfun)
cc <- lme4:::checkConv(attr(opt,"derivs"), opt$par, ctrl = lmerControl()$checkConv,
                lbound=environment(devfun)$lower)
fm1 <- mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr,
                lme4conv=cc)
fm2 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)

## basic equivalence
fm1C <- fm1
fm1C@call <- fm2@call
expect_equal(fm2,fm1C)
expect_equal(range(residuals(fm1)), c(-101.18, 132.547), tolerance = 1e-5) # these are "outliers"!
expect_is(model.frame(fm1),"data.frame")
## formulae
mfm1 <- model.frame(fm1)
expect_equal(formula(fm1),         Reaction ~ Days + (Days | Subject))
expect_equal(formula(terms(mfm1)), Reaction ~ Days + (Days + Subject))
new_form_modframe <- (getRversion() >= "3.6.0" &&
                      as.numeric(version[["svn rev"]]) >= 75891)
expect_equal(formula(mfm1),
             if(new_form_modframe) {
                 Reaction ~ Days + (Days + Subject)
             } else
                 Reaction ~ Days + Subject
             )
## predictions
expect_equal(predict(fm1,newdata=sleepstudy[1:10,],re.form=NULL),
             predict(fm2,newdata=sleepstudy[1:10,],re.form=NULL))
expect_equal(predict(fm1,newdata=sleepstudy),
             predict(fm1))

lmodOff <- lFormula(Reaction ~ Days + (Days|Subject) + offset(0.5*Days),
                    sleepstudy)
devfunOff <- do.call(mkLmerDevfun, lmodOff)
opt <- optimizeLmer(devfunOff)
fm1Off <- mkMerMod(environment(devfunOff), opt, lmodOff$reTrms, fr = lmodOff$fr)
fm2Off <- lmer(Reaction ~ Days + (Days|Subject) + offset(0.5*Days), sleepstudy)
expect_equal(predict(fm1Off,newdata=sleepstudy[1:10,],re.form=NULL),
             predict(fm2Off,newdata=sleepstudy[1:10,],re.form=NULL))

## FIXME: need more torture tests with offset specified, in different environments ...

## FIXME: drop1(.) doesn't work with modular objects ... hard to see how it
##  could, though ...
## drop1(fm1Off)
drop1(fm2Off)
