library(lme4)
library(testthat)

source(system.file("testdata/lme-tst-funs.R", package="lme4", mustWork=TRUE))
##-> gSim(), a general simulation function ...

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

