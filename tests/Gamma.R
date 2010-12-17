## generate a basic Gamma/random effects simulation
## x ~ U[0,1]
## reff ~ N(0,1)
## eta ~ 4+3*x + reff[block]
## try various links and distributions

set.seed(101)
d <- expand.grid(block=LETTERS[1:26],rep=1:100)
d$x <- runif(nrow(d))
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

## Gaussian with inverse link
dGi <- d
dGi$mu <- 1/d$eta ## inverse link
## make sd small enough to avoid negative values
dGi$y <- rnorm(nrow(d),dGi$mu,sd=0.01)

## Gamma with log link
dL <- d
dL$mu <- exp(d$eta)
dL$y <- rgamma(nrow(d),scale=dL$mu/2,shape=2)
############

library(lme4)

gm0 <- glm(y~x, data=d, family=Gamma)
gm1 <- glm(y~x+block-1, data=d, family=Gamma)
sd(coef(gm1)[-1])  ## correct
mean(coef(gm1)[-1]) ## not too bad
coef(gm1)[1]

try(gm3 <- glmer(y ~ x + (1|block), d, Gamma,
             verbose=TRUE))

## correctly specified model with "true" parameters as starting values
try(gm3 <- glmer(y ~ x + (1|block), d, Gamma,
             start=list(fixef=c(4,3),ST=list(matrix(1))),
             verbose=TRUE))
## does even worse (fails on iteration 1 with negative mu (Gauss-Newton/Fisher problem?)

###
## Poisson
gP1 <- glmer(y ~ x + (1|block), data=dP, family=poisson, verbose=TRUE)

## works just fine.

## so does Gaussian with log link
gG1 <- glmer(y ~ x + (1|block), data=dG, family=gaussian(link="log"), verbose=TRUE)

## Gaussian with inverse link
gGi1 <- glmer(y ~ x + (1|block), data=dGi, family=gaussian(link="inverse"), verbose=TRUE)
## ... sets variance to zero, converges back to GLM solution

try(gL <- glmer(y ~ x + (1|block), data=dL, family=Gamma(link="log"), verbose=TRUE))
## not so great (starts reasonable, but then goes crazy)
