library(lme4)
## generate a basic Gamma/random effects sim
set.seed(101)
d <- expand.grid(block=LETTERS[1:26],rep=1:100)
d$x <- runif(nrow(d))  ## sd=1
reff_f <- rnorm(length(levels(d$block)),sd=1)
## need intercept large enough to avoid negative values
d$eta0 <- 4+3*d$x  ## version without random effects
d$eta <- d$eta0+reff_f[d$block]
## inverse link
d$mu <- 1/d$eta
d$y <- rgamma(nrow(d),scale=d$mu/2,shape=2)

## version with Poisson errors instead (to check GLM in general)
dP <- d
dP$mu <- exp(d$eta) ## log link
dP$y <- rpois(nrow(d),dP$mu)

## version with Gaussian errors and log link (to see if there
##   might be something wrong with families that don't set scale=1?
##  need to use a non-identity link, otherwise glmer calls lmer
dG <- d
dG$mu <- exp(d$eta)
dG$y <- rnorm(nrow(d),dG$mu,sd=2)

## Gaussian with inverse link?
dGi <- d
dGi$mu <- 1/d$eta ## inverse link
## make sd small enough to avoid negative values
dGi$y <- rnorm(nrow(d),dGi$mu,sd=0.01)
############

gm0 <- glm(y~1, data=d, family=Gamma)
gm1 <- glm(y~block-1, data=d, family=Gamma)
sd(coef(gm1))

gm2 <- glmer(y ~ 1 + (1|block), d, Gamma,
             verbose=TRUE)

## do we do any better with a correctly specified model??
## no.
gm3 <- glmer(y ~ x + (1|block), d, Gamma,
             verbose=TRUE)

## correctly specified model with "true" parameters as starting values
gm3B <- glmer(y ~ x + (1|block), d, Gamma,
             start=list(fixef=c(4,3),ST=list(matrix(1))),
             verbose=TRUE)

stopifnot(all.equal(fixef(gm3),fixef(gm3B)),
          all.equal(VarCorr(gm3),VarCorr(gm3B)))

###
## Poisson
gP1 <- glmer(y ~ 1 + (1|block), data=dP, family=poisson, verbose=TRUE)
gP2 <- glmer(y ~ x + (1|block), data=dP, family=poisson, verbose=TRUE)

## works just fine.

## so does Gaussian with log link
gG1 <- glmer(y ~ 1 + (1|block), data=dG, family=gaussian(link="log"), verbose=TRUE)
gG2 <- glmer(y ~ x + (1|block), data=dG, family=gaussian(link="log"), verbose=TRUE)


## FIXME: get these working
## Gaussian with inverse link ... FIXME: not working yet
try(gGi1 <- glmer(y ~ 1 + (1|block), data=dGi, family=gaussian(link="inverse"), verbose=TRUE))
try(gGi2 <- glmer(y ~ x + (1|block), data=dGi, family=gaussian(link="inverse"), verbose=TRUE))

## sets variance to zero, converges back to GLM solution

## FIXME: cloglog link
