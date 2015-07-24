### https://github.com/lme4/lme4/issues/287
### ------------------------------------###

## I recently switched to using the development version of lme4 (1.1-8) and have run into problems with models fit using `glmer.nb` (which I know is experimental).

## I'm using an example from the `glmmadmb` help page to illustrate what I've been experiencing.

### ==> http://glmmadmb.r-forge.r-project.org/
###     =====================================
<<<<<<< Updated upstream
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
=======
install.packages("glmmADMB", repos="http://r-forge.r-project.org", type="source")

library(glmmADMB)
##      ^^^^^^^^^  not on CRAN ... arge
# Load development version of lme4
Mlibrary(lme4) ## library(lme4, lib.loc = "C:/Program Files/R/Packages")
packageVersion("lme4")
## [1] ‘1.1.8’
## Fit model with most current version of lme4 (1.1-8)
fit1 <- glmer.nb(y ~ Base*trt + Age + Visit + (1|subject), data=epil2)
## --> warning: failed to converge
fit1

## Generalized linear mixed model fit by maximum
##   likelihood (Laplace Approximation) [glmerMod]
##  Family: Negative Binomial(0.65)  ( log )
## Formula: y ~ Base * trt + Age + Visit + (1 | subject)
##    Data: ..2
##       AIC       BIC    logLik  deviance  df.resid
## 1421.5207 1449.2314 -702.7604 1405.5207       228
## Random effects:
##  Groups   Name        Std.Dev.
##  subject  (Intercept) 1.171e-06
##  Residual             6.109e-01
## Number of obs: 236, groups:  subject, 59
## Fixed Effects:
##       (Intercept)               Base       trtprogabide
##           -1.2401             0.8790            -0.8726
##               Age              Visit  Base:trtprogabide
##            0.4762            -0.2725             0.3342


## The estimate of theta is very small and the random effects variance is essentially 0.
## Compare this to the output from a model fit with `glmmadmb`.


## Compare with glmmadmb fit
epil2$subject <-  factor(epil2$subject)
fit2 <- glmmadmb(y ~ Base*trt + Age + Visit + (1|subject), data=epil2, family="nbinom")
fit2
>>>>>>> Stashed changes

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

<<<<<<< Updated upstream
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

=======
## Things look more reasonable to me if I go back to lme4 version 1.1-7.
##--------------------------------------------------------------------------

## Detach package lme4 version 1.1-8 and load 1.1-7
detach("package:lme4", unload=TRUE)
library(lme4, lib.loc = "~/R/Pkgs/lme4.Rcheck-64b.bakk" ## <- nb-mm4
        )## "C:/Program Files/R/R-3.1.3/library")
## check:
stopifnot(packageVersion("lme4") == "1.1.7")

# Refit model with lme4 1.1-7
fit3 <- glmer.nb(y ~ Base*trt + Age + Visit + (1|subject), data=epil2)
fit3

vc3 <- VarCorr(fit3)
stopifnot(all.equal(c(as.vector(vc3$subject), attr(vc3, "sc")),
                    c(0.202958669916721, 0.968798378253818), tol=4e-15))

dput(fixef(fit3))
fixef.Exp <- c("(Intercept)" = -1.33798122559697,
               Base = 0.882641577389333, trtprogabide = -0.927943548344718,
               Age = 0.474584264700356, Visit = -0.268534514764688,
               `Base:trtprogabide` = 0.336536718825345)
stopifnot(all.equal(fixef(fit3), fixef.Exp, tol=1e-14))

## Generalized linear mixed model fit by maximum
##   likelihood (Laplace Approximation) [glmerMod]
##  Family: Negative Binomial(7.4734)  ( log )
## Formula: y ~ Base * trt + Age + Visit + (1 | subject)
##    Data: ..2
##       AIC       BIC    logLik  deviance  df.resid
## 1264.9701 1292.6808 -624.4851 1248.9701       228
## Random effects:
##  Groups   Name        Std.Dev.
##  subject  (Intercept) 0.4505
##  Residual             0.9688
## Number of obs: 236, groups:  subject, 59
## Fixed Effects:
##       (Intercept)               Base       trtprogabide
##           -1.3380             0.8826            -0.9279
##               Age              Visit  Base:trtprogabide
##            0.4746            -0.2685             0.3365

## After poking around for a bit, my best guess is that the `deviance.merMod`
## function in the development version of lme4
## now returns the sum of the squared deviance residuals for GLMM's,
## where it used to return twice the negative log-likelihood.
## The `deviance` function is used in the final theta optimization step of `glmer.nb`
## (in the `optTheta` function).
>>>>>>> Stashed changes
