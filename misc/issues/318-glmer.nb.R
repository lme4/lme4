## From: Ben Bolker <notifications@github.com>
## Subject: [lme4] interaction between glmer.nb and effects package (#318)
## Date: Thu, 23 Jul 2015 07:28:14 -0700

## Let's see if we can hack this on our side by putting `th` somewhere in the environment: if not we may need to ask the `effects` maintainers for help.


set.seed(101)
dd <- expand.grid(f1=factor(1:3),f2=factor(1:3),g=1:10,rep=1:10)
dd$y <- rnbinom(nrow(dd), mu=2, size=0.5)
library(lme4)
m.nb <- glmer.nb(y ~ f1*f2 + (1|g), data=dd, verbose = 2) ## >> *TONS* of output
##                                         ^^^^^^^^^^^^ new
## MM: here (inside glmer.nb !), we need to replace 'th' by the number
##    1) that's done.
## ==> need to replace '..2'  by 'dd' ... which is a bit more work

## Issue #266 -- already has same problem:
plot(m.nb, g ~ resid(.))

## _AND_  issue #176 :
m2 <- update(m.nb)
## but the error here is more strange: 'verbose' ?????
## Error in lme4::glFormula(formula = y ~ f1 * f2 + (1 | g), data = ..2,  :
##   object 'verbose' not found
##
## specifying 'verbose' explitly helps
##
m2 <- update(m.nb, verbose=FALSE)
## Error: 'data' not found, and some variables missing from formula environment
## ... so it may just need to get the real 'data', indeed,
##     it works  if you specify 'data' yourself :
m2 <- update(m.nb, data=dd, verbose=FALSE)



library(effects)
library(MASS)
effect("f1*f2",m.nb)

## After the 'th' fix :
## Error in is.data.frame(data) :
##   ..2 used in an incorrect context, no ... to look in

traceback()
## 13: is.data.frame(data)
## 12: model.frame.default(formula = y ~ f1 * f2, data = ..2, drop.unused.levels = TRUE)
## 11: stats::model.frame(formula = y ~ f1 * f2, data = ..2, drop.unused.levels = TRUE)
## 10: eval(expr, envir, enclos)
## 9: eval(mf, parent.frame())
## 8: glm(formula = y ~ f1 * f2, family = negative.binomial(theta = 0.0244908267262883),
##        data = ..2)
##_______________^^^___________________
## 7: eval(expr, envir, enclos)
## 6: eval(cl)
## 5: mer.to.glm(mod, KR = KR)
## 4: effect(term, mer.to.glm(mod, KR = KR), vcov., ...)
## 3: effect.mer(term, mod, vcov. = vcov, KR = KR, ...)
## 2: effect.merMod("f1*f2", m.nb)
## 1: effect("f1*f2", m.nb)




## The above reveals a small glitch
## in the output of MASS::theta.ml(*, trace=TRUE):
require(MASS)
quine.nb <- glm.nb(Days ~ .^2, data = quine)
##
theta.ml(quine$Days, fitted(quine.nb), trace=TRUE)
##
## theta.ml: iter 0 'theta = 1.574510'  <-- Has extra quotes
## theta.ml: iter1 theta =1.60293       <-- space missing after 'iter'
## theta.ml: iter2 theta =1.60364
## theta.ml: iter3 theta =1.60364

theta.ml(quine$Days, fitted(quine.nb), trace=TRUE, eps = 1e-10)
## Here, the user might want to see slightly more digits
