library("lme4")

## tests of a variety of GLMM families and links
## coding: family {g=Gamma, P=Poisson, G=Gaussian, B=binomial}
##         link   {l=log, i=inverse, c=cloglog, i=identity}
##         model  {1 = intercept-only, 2 = with continuous predictor}

set.seed(101)
d <- expand.grid(block=LETTERS[1:26], rep=1:100, KEEP.OUT.ATTRS = FALSE)
d$x <- runif(nrow(d))  ## sd=1
reff_f <- rnorm(length(levels(d$block)),sd=1)
## need intercept large enough to avoid negative values
d$eta0 <- 4+3*d$x  ## fixed effects only
d$eta <- d$eta0+reff_f[d$block]

## Gamma, inverse link
d$mu <- 1/d$eta
d$y <- rgamma(nrow(d),scale=d$mu/2,shape=2)

## Gamma, log link
dgl <- d
dgl$mu <- exp(d$eta)
dgl$y <- rgamma(nrow(d),scale=dgl$mu/2,shape=2)

## Poisson, log link
dP <- d
dP$mu <- exp(d$eta) ## log link
dP$y <- rpois(nrow(d),dP$mu)

##  Gaussian, log link
##  need to use a non-identity link, otherwise glmer calls lmer
dG <- d
dG$mu <- exp(d$eta)
dG$y <- rnorm(nrow(d),dG$mu,sd=2)

## Gaussian with inverse link
dGi <- d
dGi$mu <- 1/d$eta ## inverse link
## make sd small enough to avoid negative values
dGi$y <- rnorm(nrow(d),dGi$mu,sd=0.01)

## binomial with cloglog link
dBc <- d
cc <- binomial(link="cloglog")
dBc$mu <- cc$linkinv(d$eta - 5)         # -5, otherwise y will be constant
dBc$y <- factor(rbinom(nrow(d),dBc$mu,size=1))

## binomial with identity link
dBi <- d
cc <- binomial(link="identity")
dBi$mu <- cc$linkinv(d$eta/10)         # scale so range goes from 0.2-0.8
dBi$y <- factor(rbinom(nrow(d),dBi$mu,size=1))


############
## Gamma/inverse

## GLMs
gm0 <- glm(y ~ 1,       data=d, family=Gamma)
gm1 <- glm(y ~ block-1, data=d, family=Gamma)
stopifnot(all.equal(sd(coef(gm1)),1.00753942148611))

gm2 <- glmer(y ~ 1 + (1|block), d, Gamma, nAGQ=0)
gm3 <- glmer(y ~ x + (1|block), d, Gamma, nAGQ=0)
gm2B <- glmer(y ~ 1 + (1|block), d, Gamma)
gm3B <- glmer(y ~ x + (1|block), d, Gamma)

##
## library(hglm)
## h1 <- hglm2(y~x+(1|block), data=d, family=Gamma())
## lme4.0 fails on all of these ...

## Gamma/log
ggl1 <- glmer(y ~ 1 + (1|block), data=dgl, family=Gamma(link="log"))
ggl2 <- glmer(y ~ x + (1|block), data=dgl, family=Gamma(link="log"))

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
gGi1 <- glmer(y ~ 1 + (1|block), data=dGi,family=gaussian(link="inverse"))
gGi2 <- glmer(y ~ x + (1|block), data=dGi, family=gaussian(link="inverse"))


## Binomial/cloglog
gBc1 <- glmer(y ~ 1 + (1|block), data=dBc, family=binomial(link="cloglog"))

gBc2 <- glmer(y ~ x + (1|block), data=dBc,
              family=binomial(link="cloglog"))
## library("glmmADMB")
## glmmadmbfit <- glmmadmb(y ~ x + (1|block), data=dBc,
## family="binomial",link="cloglog")
glmmadmbfit <- structure(list(fixef = structure(c(-0.717146132730349, 2.83642900561633),
                        .Names = c("(Intercept)", "x")), VarCorr = structure(list(
                        block = structure(0.79992, .Dim = c(1L, 1L),
                        .Dimnames = list(
                                   "(Intercept)", "(Intercept)"))),
                           .Names = "block", class = "VarCorr")),
                            .Names = c("fixef", "VarCorr"))
stopifnot(all.equal(fixef(gBc2),glmmadmbfit$fixef,tolerance=5e-3))
## pretty loose tolerance ...
stopifnot(all.equal(unname(unlist(VarCorr(gBc2))),
                    c(glmmadmbfit$VarCorr$block),tolerance=2e-2))

gBi1 <- glmer(y ~ 1 + (1|block), data=dBi, family=binomial(link="identity"))
gBi2 <- glmer(y ~ x + (1|block), data=dBi, family=binomial(link="identity"))

## FIXME: should test more of the *results* of these efforts, not
##  just that they run without crashing ...
