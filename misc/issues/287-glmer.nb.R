### https://github.com/lme4/lme4/issues/287
### ------------------------------------###

## I recently switched to using the development version of lme4 (1.1-8) and have run into problems with models fit using `glmer.nb` (which I know is experimental).

## I'm using an example from the `glmmadmb` help page to illustrate what I've been experiencing.

### ==> http://glmmadmb.r-forge.r-project.org/
###     =====================================
if(!require("glmmADMB"))
    install.packages("glmmADMB", repos="http://r-forge.r-project.org", type="source")

library(glmmADMB)
##      ^^^^^^^^^  not on CRAN, (but R-forge, see above)
packageVersion("glmmADMB")

# Load development version of lme4
library(lme4) ## library(lme4, lib.loc = "C:/Program Files/R/Packages")
packageVersion("lme4")
## [1] ‘1.1.9’ (<- now)
fit4 <- glmer.nb(y ~ Base*trt + Age + Visit + (1|subject), data = epil2)
## --> warnings: failed to converge

## Compare with glmmadmb fit
epil2$subject <-  factor(epil2$subject)
fitAD <- glmmadmb(y ~ Base*trt + Age + Visit + (1|subject), data = epil2,
                  family="nbinom")
fitAD

## GLMM's in R powered by AD Model Builder:

##   Family: nbinom
##   alpha = 7.4702
##   link = log

## Fixed effects:
##   Log-likelihood: -624.551
##   AIC: 1265.102
##   Formula: y ~ Base * trt + Age + Visit + (1 | subject)
##       (Intercept)              Base      trtprogabide
##        -1.3300404         0.8839167        -0.9299655
##               Age             Visit Base:trtprogabide
##         0.4751431        -0.2701602         0.3372419

## Random effects:
## Structure: Diagonal matrix
## Group=subject
##             Variance StdDev
## (Intercept)   0.2172  0.466

## Number of observations: total=236, subject=59
##--------------------------------------------------------------------------

all.equal(fixef(fit4), coef(fitAD))
## [1] "Mean relative difference: 0.003201962"
stopifnot(all.equal(fixef(fit4),
                    coef (fitAD), tol = 0.008))
VarCorr(fit4)
## Groups   Name       Std.Dev.
## subject (Intercept) 0.46579

VarCorr(fitAD)
## Group=subject
##             Variance StdDev
## (Intercept)   0.2172  0.466

all.equal(VarCorr(fit4 )[["subject"]][1,1],
          VarCorr(fitAD)[["subject"]][1,1]) # 0.000955
stopifnot(all.equal(VarCorr(fit4 )[["subject"]][1,1],
                    VarCorr(fitAD)[["subject"]][1,1], tol= 0.002))


all.equal(data.matrix(ranef(fit4 )[["subject"]]),
                      ranef(fitAD)[["subject"]])# 0.028575..
stopifnot(all.equal(data.matrix(ranef(fit4 )[["subject"]]),
                                ranef(fitAD)[["subject"]],
                    tol = 0.05))

## The negative binomial shape parameter 'theta'
cat(sprintf("%.10g (%.10g)\n", fitAD$alpha, fitAD$sd_alpha))
## -> 7.4702 (1.7605)
(th.4 <- getME(fit4, "glmer.nb.theta")) # 7.456566

all.equal(th.4, fitAD$alpha) # 0.001828
stopifnot(all.equal(th.4, fitAD$alpha, tol = 0.004))

