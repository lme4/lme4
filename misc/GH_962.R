form <- out_bin_1 ~ ftime + trt + ar1(ftime + 0 | id_cluster) +
  (1 | id_individual)

form2 <- update(form, . ~ . - ar1(ftime + 0 | id_cluster) +
                        cs(ftime + 0 | id_cluster))

form3 <- update(form, . ~ . - ar1(ftime + 0 | id_cluster) +
                        (1|id_cluster/ftime))

form4 <- update(form, . ~ . - (1 | id_individual))

library(glmmTMB)
library(lme4)
stopifnot(packageVersion("lme4") >= "2.0.2")
set.seed(101)
sfun <- function(n) factor(sample(1:n, replace=TRUE, size = 1465))
dd <- data.frame(id_cluster = sfun(7), id_individual = sfun(486),
                 ftime = sfun(4), trt = sfun(2))
np <- list(beta = c(0, rep(0.1, 3), 1), theta = rep(0, 3))
dd$out_bin_1 <- simulate_new(form[-2], newdata = dd,
             newparams = np,
             family = binomial)[[1]]

## devtools::load_all()
## debug(optimizeGlmer)
## looks OK after nAGQ0: environment(devfun)$pp$beta(1)
## but not after nAGQ1

## looks like it's within optimizeGlmer -> devfun() -> pwrssUpdate() ?
## inside .Call(glmerLaplace, ...) ?

## + + *** pwrssUpdate step 0
## pdev=1351.8; delu_min: -1.52241; delu_max: 1.69336; delb_min: 0; delb_max: 0
## *** pwrssUpdate step 1
## pdev=1351.8; delu_min: -1.52241; delu_max: 1.69336; delb_min: 0; delb_max: 0
## [1] 1638.19

## ?? set debug flag internally?
## try git bisect ... ??
## narrow down cases where this happens? (is it true for *all* 

## Not universal. e.g.

gm.cs <- glmer(use ~ age + urban + cs(1 + urban | district),
               data = Contraception,
               family = binomial)
fixef(gm.cs)  ## OK

gm.ar1_0 <- glmer(use ~ age + urban + ar1(0 + urban | district),
               data = Contraception,
               family = binomial)
colnames(getME(gm.ar1_0, "X"))
getME(gm.ar1_0, "beta")  ## only two values ???

try(fixef(gm.ar1_0))
## Error in attributes(.Data) <- c(attributes(.Data), attrib) (from lmer.R#971) : 
##   'names' attribute [3] must be the same length as the vector [2]



## but ...
gm.ar1 <- glmer(use ~ urban + ar1(0 + factor(age) | district),
               data = Contraception,
               family = binomial)
## uh-oh
fixef(gm.ar1)

## still not sure what's going on ...
## cross-check with other forms?

## from a833370f0b10ad221b43c400fa1 (pre-2.0), for comparison
## devfun for glmer/nAGQ > 0
olddevfun <- function(pars) {
  ## pp$setDelu(rep(0, length(pp$delu)))
  resp$setOffset(baseOffset)
  resp$updateMu(lp0)
  pp$setTheta(as.double(pars[dpars])) # theta is first part of pars
  spars <- as.numeric(pars[-dpars])
  offset <- if (length(spars)==0) baseOffset else baseOffset + pp$X %*% spars
  resp$setOffset(offset)
  p <- pwrssUpdate(pp, resp, tol=tolPwrss, GQmat=GQmat,
                   compDev=compDev, grpFac=fac, maxit=maxit, verbose=verbose)
  resp$updateWts()
  p
}

lme4_1 <- glmer(form, data = dd, family = binomial)
fixef(lme4_1)
lme4_1@optinfo$val

lme4_1B <- update(lme4_1, nAGQ=0)
fixef(lme4_1B)

lme4_2 <- update(lme4_1, form = form2)

lme4_2@optinfo$val

lme4_3 <- update(lme4_1, form = form3)
fixef(lme4_3)
lme4_3@optinfo$val

## ar1 only, no individual RE
lme4_4 <- update(lme4_1, form = form4)
fixef(lme4_4)

## only happens for 'new-style' (structured); something happens
## (or fails to happen) in the parameter mapping step?

glmmTMB_1 <- glmmTMB(form, data = dd, family = binomial)
c(logLik(lme4_1) - logLik(glmmTMB_1))
logLik(lme4_1)
logLik(glmmTMB_1)


dd$out_gauss <- simulate_new(form[-2], newdata = dd,
             newparams = c(np, list(betadisp=0)),
             family = gaussian)[[1]]

form2 <- update(form, out_gauss ~ .)
lme4_2 <- lmer(form2, data = dd, REML = FALSE)
glmmTMB_2 <- glmmTMB(form2, data = dd, family = gaussian)
fixef(lme4_2)
