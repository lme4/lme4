library(lme4)
library(testthat)

source(system.file("testdata/lme-tst-funs.R", package="lme4", mustWork=TRUE))
##-> gSim(), a general simulation function ...

## hand-coded Pearson residuals {for  sumFun() }
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
d.P <- gSim(family=poisson())
g2 <-  glmer(y ~ x + (1|block),data = d.P, family=poisson)
expect_equal(var(sumFun(g2)),0)

g2W <- glmer(y ~ x + (1|block), data=d.P,
             family=poisson, weights=rep(2,nrow(d.P)))
expect_equal(var(sumFun(g2W)),0)
## correct

## non-Bernoulli BINOMIAL
g3 <- glmer(y ~ x + (1|block),data=gSim(family=binomial(), nbinom=10),
            family=binomial)
expect_equal(var(sumFun(g3)),0)

d.b.2 <- gSim(nperblk = 2, family=binomial())
g.b.2 <- glmer(y ~ x + (1|block), data=d.b.2, family=binomial)
expect_equal(var(sumFun(g.b.2)), 0)

## Many blocks of only 2 observations each - (but nicely balanced)
## Want this "as" https://github.com/lme4/lme4/issues/47
## (but it "FAILS" survival already):
##
## n2 = n/2 :
n2 <- 2048
n2 <-  100 # for building/testing
set.seed(47)
dB2 <- gSim(n2, nperblk = 2, x= rep(0:1, each= n2), family=binomial())
##                       --  --     ---  --------
gB2 <- glmer(y ~ x + (1|block), data=dB2, family=binomial)

## FAILS -----
## library(survival)
## (gSurv.B2 <- clogit(y ~ x + strata(block), data=dB2))
## summary(gSurv.B2)
## (SE.surf <- sqrt(diag(vcov(gSurv.B2))))



g3 <-  glmer(y ~ x + (1|block),data=gSim(family=binomial(),nbinom=10),
             family=binomial)
expect_equal(var(sumFun(g3)),0)

## check dispersion parameter
## (lowered tolerance to pass checks on my machine -- SCW)
expect_equal(sigma(g0)^2,0.4888248,tol=1e-4)

