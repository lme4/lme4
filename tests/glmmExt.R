library(lme4)

## tests of a variety of GLMM families and links

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
str(d)## 2600 obs .. 'block' with 26 levels

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

############
## Gamma/inverse

## GLMs
gm0 <- glm(y ~ 1,       data=d, family=Gamma)
gm1 <- glm(y ~ block-1, data=d, family=Gamma)
sd(coef(gm1)) # 1.007539

gm2 <- glmer(y ~ 1 + (1|block), d, Gamma, verbose = 4)
gm3 <- glmer(y ~ x + (1|block), d, Gamma, verbose = 4)

## with "true" parameters as starting values
gm3B <- glmer(y ~ x + (1|block), d, Gamma,
             start=list(fixef=c(4,3),ST=list(matrix(1))),
             verbose = 4)

stopifnot(all.equal(fixef  (gm3),fixef  (gm3B)),
          all.equal(VarCorr(gm3),VarCorr(gm3B)))

## Gamma/log
ggl1 <- glmer(y ~ 1 + (1|block), data=dgl, family=Gamma(link="log"), verbose= 2)
ggl2 <- glmer(y ~ x + (1|block), data=dgl, family=Gamma(link="log"), verbose= 2)

##
## library(lme4.0)
## ggl1 <- glmer(y ~ 1 + (1|block), data=dgl, family=Gamma(link="log"), verbose= 2)
## fails

## Poisson/log
gP1 <- glmer(y ~ 1 + (1|block), data=dP, family=poisson, verbose= 2)
gP2 <- glmer(y ~ x + (1|block), data=dP, family=poisson, verbose= 2)

## Gaussian/log
gG1 <- glmer(y ~ 1 + (1|block), data=dG, family=gaussian(link="log"), verbose=TRUE)
gG2 <- glmer(y ~ x + (1|block), data=dG, family=gaussian(link="log"), verbose=TRUE)

## works with lme4.0 but AIC/BIC/logLik are crazy, and scale
## parameter is not reported
## glmmML etc. doesn't allow models with scale parameters
## gG1B <- glmmadmb(y ~ 1 + (1|block), data=dG,
##                  family="gaussian",link="log",verbose=TRUE)
## what is the best guess at the estimate of the scale parameter?
## is it the same as sigma?
## gG1B$alpha

if(Sys.info()["user"] != "maechler") { # <- seg.faults (MM)

## FIXME: fail for MM (retest?)
## Gaussian/inverse
    gGi1 <- glmer(y ~ 1 + (1|block), data=dGi, family=gaussian(link="inverse"), verbose= 3)
    gGi2 <- glmer(y ~ x + (1|block), data=dGi, family=gaussian(link="inverse"), verbose= 3)
}


## Binomial/cloglog
gBc1 <- glmer(y ~ 1 + (1|block), data=dBc,
              family=binomial(link="cloglog"), verbose= 3)
if (FALSE)                              # still having problems with this one
    gBc2 <- glmer(y ~ x + (1|block), data=dBc,
                  family=binomial(link="cloglog"), verbose= 3)
