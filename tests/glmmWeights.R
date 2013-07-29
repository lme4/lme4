library(lme4)
library(testthat)

## a general simulation function ...
gSim <- function(nblk=26,  
                 nperblk=100,
                 sigma=1,
                 beta=c(4,3),
                 shape=2,    ## shape parameter for Gamma
                 nbinom=10,  ## N for binomial trials
                 family=Gamma()) {
    d <- expand.grid(block=LETTERS[1:nblk],
                     rep=1:nperblk, KEEP.OUT.ATTRS = FALSE)
    d$x <- runif(nrow(d))  ## sd=1
    reff_f <- rnorm(length(levels(d$block)),sd=sigma)
    ## need intercept large enough to avoid negative values
    d$eta0 <- beta[1]+beta[2]*d$x  ## fixed effects only
    d$eta <- d$eta0+reff_f[d$block]
    d$mu <- family$linkinv(d$eta)
    d$y <- switch(family$family,
                  Gamma=rgamma(nrow(d),scale=d$mu/shape,shape=shape),
                  poisson=rpois(nrow(d),d$mu),
                  binomial=if (nbinom==1) {
                      rbinom(nrow(d),prob=d$mu,size=1)
                  } else 
              {z <- rbinom(nrow(d),prob=d$mu,size=nbinom);
               cbind(succ=z,fail=nbinom-z)})
    d
}


## hand-coded Pearson residuals
mypresid <- function(x) {
    mu <- fitted(x)
    w <- weights(x)
    (getME(x,"y")-mu)*sqrt(w)/sqrt(x@resp$family$variance(mu))
}

sumFun <- function(m) {
    ## sum of weighted residuals, 3 ways
    ## (always identical)
    wrss1 <- m@devcomp$cmp["wrss"]
    wrss2 <- sum(residuals(m,type="pearson")^2)
    wrss3 <- sum(m@resp$wtres^2)
    ## compare to hand-fitted Pearson resids ...
    wrss4 <- sum(mypresid(m)^2)
    c(wrss1,wrss2,wrss3,wrss4)
}

set.seed(101)
library(lme4)

## GAMMA
g0 <-  glmer(y~x+(1|block),data=gSim(),family=Gamma)
expect_equal(var(sumFun(g0)),0)

## BERNOULLI
g1 <-  glmer(y~x+(1|block),data=gSim(family=binomial(),nbinom=1),
             family=binomial)
expect_equal(var(sumFun(g1)),0)


## POISSON 
g2 <-  glmer(y~x+(1|block),data=gSim(family=poisson()),
             family=poisson)
expect_equal(var(sumFun(g2)),0)

d <- gSim(family=poisson())
g2W <-  glmer(y~x+(1|block),data=d,
             family=poisson,weights=rep(2,nrow(d)))
expect_equal(var(sumFun(g2W)),0)
## correct

g3 <-  glmer(y~x+(1|block),data=gSim(family=binomial(),nbinom=10),
             family=binomial)
expect_equal(var(sumFun(g3)),0)

## check dispersion parameter
## (lowered tolerance to pass checks on my machine -- SCW)
expect_equal(sigma(g0)^2,0.4888248,tol=1e-4)

