## Tests of a variety of GLMM families and links
## coding: family {g=Gamma, P=Poisson, G=Gaussian, B=binomial}
##         link   {l=log, i=inverse, c=cloglog, i=identity}
##         model  {1 = intercept-only, 2 = with continuous predictor}

library("lme4")

source(system.file("testdata/lme-tst-funs.R", package="lme4", mustWork=TRUE))
##-> gSim(), a general simulation function ...

str(gSim)
## function (nblk = 26, nperblk = 100, sigma = 1, beta = c(4, 3),
##           x = runif(n), shape = 2, nbinom = 10, family = Gamma())

if (.Platform$OS.type != "windows") withAutoprint({
set.seed(101)
## Gamma, inverse link (= default) :
d <- gSim()
## Gamma, log link   eta = log(mu) :
dgl <- gSim(dInitial = d, family = Gamma(link = log))
## Poisson, log link
dP <- gSim(dInitial = d, family = poisson())
##  Gaussian, log link --- need to use a non-identity link, otherwise glmer calls lmer
dG <- gSim(dInitial = d, family = gaussian(link = log), sd = 2)
## Gaussian with inverse link :     (sd small enough to avoid negative values) :
dGi <- gSim(dInitial = d, family = gaussian(link = inverse), sd = 0.01)
## binomial with cloglog link
dBc <- d
dBc$eta <- d$eta - 5 # <==> beta intercept 5 less: otherwise  y will be constant
dBc <- gSim(dInitial = dBc, ## beta = c(-1, 3),
            nbinom = 1, family = binomial(link="cloglog"))

## binomial with identity link
dBi <- d
dBc$eta <- d$eta / 10 # <==> beta slope / 10 : scale so range goes from 0.2-0.8
dBi <- gSim(dInitial = dBc, ## beta = c(4, 3/10),
            nbinom = 1, family = binomial(link="identity"))


############
## Gamma/inverse

## GLMs
gm0 <- glm(y ~ 1,       data=d, family=Gamma)
gm1 <- glm(y ~ block-1, data=d, family=Gamma)
stopifnot(all.equal(sd(coef(gm1)),1.00753942148611))

gm2  <- glmer(y ~ 1 + (1|block), d, Gamma, nAGQ=0)
gm3  <- glmer(y ~ x + (1|block), d, Gamma, nAGQ=0)
gm2B <- glmer(y ~ 1 + (1|block), d, Gamma)
gm3B <- glmer(y ~ x + (1|block), d, Gamma)

##   y ~ x + (1|block),  Gamma   is TRUE model
summary(gm3)
summary(gm3B)# should be better
## Both have "correct" beta ~= (4, 3)  -- but *too* small  (sigma_B, sigma) !!
## reference values updated for the dispersion-profiling fix (phi is now
## profiled rather than fixed at 1 for free-dispersion families like
## Gamma), and again for disp_dof_correction=TRUE (now the default),
## which applies an n-qEff degrees-of-freedom correction on top -- see
## misc/Gamma_GLMM/README_Gamma_GLMMs.md: fixed effects are largely, but
## not exactly, unaffected by either fix (both perturb the joint
## estimation slightly, via PIRLS's working weights)
stopifnot(exprs = {
    all.equal(fixef(gm3 ), c(`(Intercept)` = 4.0923502324289, x = 3.07339575816978), tol = 1e-5)
    all.equal(fixef(gm3B), c(`(Intercept)` = 4.12612735957073, x = 3.05914408107665), tol = 1e-5)
})
VarCorr(gm3)  # both variances / std.dev. should be ~ 1  but are too small



##
## library(hglm)
## h1 <- hglm2(y~x+(1|block), data=d, family=Gamma())
## lme4.0 fails on all of these ...

## Gamma/log
ggl1 <- glmer(y ~ 1 + (1|block), data=dgl, family=Gamma(link="log"))
ggl2 <- glmer(y ~ x + (1|block), data=dgl, family=Gamma(link="log"))# true model
(h.1.2 <- anova(ggl1, ggl2))
## reference values updated for the dispersion-profiling fix (phi is now
## profiled rather than fixed at 1 for free-dispersion families like
## Gamma), and again for disp_dof_correction=TRUE (now the default),
## which applies an n-qEff degrees-of-freedom correction on top -- see
## misc/Gamma_GLMM/README_Gamma_GLMMs.md
stopifnot(
    all.equal(unlist(h.1.2[2,]),
              c(npar = 4, AIC = 34230.0757742238, BIC = 34253.5288411198,
                logLik = -17111.0378871119, "-2*log(L)" = 34222.0757742238,
                Chisq = 2437.21070033035, Df = 1, `Pr(>Chisq)` = 0))
)
## "true" model :
summary(ggl2)
VarCorr(ggl2)

##
## library(lme4.0)
## ggl1 <- glmer(y ~ 1 + (1|block), data=dgl, family=Gamma(link="log"), verbose= 2)
## fails

## Poisson/log
gP1 <- glmer(y ~ 1 + (1|block), data=dP, family=poisson)
gP2 <- glmer(y ~ x + (1|block), data=dP, family=poisson)

## Gaussian/log
gG1 <- glmer(y ~ 1 + (1|block), data=dG, family=gaussian(link="log"))
gG2 <- glmer(y ~ x + (1|block), data=dG, family=gaussian(link="log"))

## works with lme4.0 but AIC/BIC/logLik are crazy, and scale
## parameter is not reported
## glmmML etc. doesn't allow models with scale parameters
## gG1B <- glmmadmb(y ~ 1 + (1|block), data=dG,
##                  family="gaussian",link="log",verbose=TRUE)
## what is the best guess at the estimate of the scale parameter?
## is it the same as sigma?
## gG1B$alpha

## if(Sys.info()["user"] != "maechler") { # <- seg.faults (MM)

## Gaussian/inverse
gGi1 <- glmer(y ~ 1 + (1|block), data=dGi, family=gaussian(link="inverse"))
gGi2 <- glmer(y ~ x + (1|block), data=dGi, family=gaussian(link="inverse"))


## Binomial/cloglog
gBc1 <- glmer(y ~ 1 + (1|block), data=dBc, family=binomial(link="cloglog"))
gBc2 <- glmer(y ~ x + (1|block), data=dBc, family=binomial(link="cloglog"))
## library("glmmADMB")
## glmmadmbfit <- glmmadmb(y ~ x + (1|block), data=dBc,
## family="binomial",link="cloglog")
glmmadmbfit <- list(fixef = c("(Intercept)" = -0.717146132730349, x =2.83642900561633),
                    VarCorr = structure(list(
                        block = structure(0.79992, dim = c(1L, 1L),
                                          dimnames = list("(Intercept)", "(Intercept)"))),
                        class = "VarCorr"))
stopifnot(all.equal(fixef(gBc2), glmmadmbfit$fixef, tolerance=5e-3))
## pretty loose tolerance ...
stopifnot(all.equal(unname(unlist(VarCorr(gBc2))),
                    c(glmmadmbfit$VarCorr$block), tolerance=2e-2))

gBi1 <- glmer(y ~ 1 + (1|block), data=dBi, family=binomial(link="identity"))
gBi2 <- glmer(y ~ x + (1|block), data=dBi, family=binomial(link="identity"))

## FIXME: should test more of the *results* of these efforts, not
##  just that they run without crashing ...
}) ## skip on windows (for speed)
