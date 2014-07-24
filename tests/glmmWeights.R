library(lme4)
library(testthat)

source(system.file("testdata/lme-tst-funs.R", package="lme4", mustWork=TRUE))
##-> gSim(), a general simulation function ...

## hand-coded Pearson residuals {for  sumFun() }
mypresid <- function(x) {
    mu <- fitted(x)
    (getME(x,"y") - mu) * sqrt(weights(x)) / sqrt(x@resp$family$variance(mu))
}

## should be equal (up to numerical error) to weights(.,type="working")
workingWeights <- function(mod) mod@resp$weights*(mod@resp$muEta()^2)/mod@resp$variance()

##' Sum of weighted residuals, 4 ways; the last three are identical
sumFun <- function(m) {
    wrss1 <- m@devcomp$cmp["wrss"]
    wrss2 <- sum(residuals(m,type="pearson")^2)
    wrss3 <- sum(m@resp$wtres^2)
    ## compare to hand-fitted Pearson resids ...
    wrss4 <- sum(mypresid(m)^2)
    c(wrss1,wrss2,wrss3,wrss4)
}
## The relative "error"/differences of the weights w[] entries
rel.diff <- function(w) abs(1 - w[-1]/w[1])

set.seed(101)

## GAMMA
g0 <-  glmer(y~x+(1|block),data=gSim(),family=Gamma)
expect_true(all(rel.diff(sumFun(g0)) < 1e-13))
expect_equal(weights(g0, type = "working"), workingWeights(g0),
             tolerance = 1e-4)  ## FIXME:  why is such a high tolerance required?

## BERNOULLI
g1 <-  glmer(y~x+(1|block),data=gSim(family=binomial(),nbinom=1),
             family=binomial)
expect_true(all(rel.diff(sumFun(g1)) < 1e-13))
expect_equal(weights(g1, type = "working"), workingWeights(g1),
             tolerance = 1e-5)  ## FIXME:  why is such a high tolerance required?


## POISSON
(n <- nrow(d.P <- gSim(family=poisson())))
g2 <-  glmer(y ~ x + (1|block), data = d.P, family=poisson)
g2W <- glmer(y ~ x + (1|block), data = d.P, family=poisson, weights = rep(2,n))
expect_true(all(rel.diff(sumFun(g2 )) < 1e-13))
expect_true(all(rel.diff(sumFun(g2W)) < 1e-13))
## correct
expect_equal(weights(g2, type = "working"), workingWeights(g2),
             tolerance = 1e-5)  ## FIXME:  why is such a high tolerance required?
expect_equal(weights(g2W, type = "working"), workingWeights(g2W),
             tolerance = 1e-5)  ## FIXME:  why is such a high tolerance required?


## non-Bernoulli BINOMIAL
g3 <- glmer(y ~ x + (1|block), data= gSim(family=binomial(), nbinom=10),
            family=binomial)
expect_true(all(rel.diff(sumFun(g3)) < 1e-13))
expect_equal(weights(g3, type = "working"), workingWeights(g3),
             tolerance = 1e-4)  ## FIXME:  why is such a high tolerance required?



d.b.2 <- gSim(nperblk = 2, family=binomial())
g.b.2 <- glmer(y ~ x + (1|block), data=d.b.2, family=binomial)

expect_true(all(rel.diff(sumFun(g.b.2 )) < 1e-13))


## Many blocks of only 2 observations each - (but nicely balanced)
## Want this "as" https://github.com/lme4/lme4/issues/47
## (but it "FAILS" survival already):
##
## n2 = n/2 :
n2 <- 2048
if(FALSE)
n2 <-  100 # for building/testing
set.seed(47)
dB2 <- gSim(n2, nperblk = 2, x= rep(0:1, each= n2), family=binomial())
##                       --  --     ---  --------
gB2 <- glmer(y ~ x + (1|block), data=dB2, family=binomial)
expect_true(all(rel.diff(sumFun(gB2)) < 1e-13))

## NB: Finite sample bias of  \hat\sigma_1 and  \hat\beta_1 ("Intercept")
##     tend to zero only slowly for  n2 -> Inf,  e.g., for
## n2 = 2048,  b1 ~= 4.3 (instead of 4);  s1 ~= 1.3 (instead of 1)

## FAILS -----
## library(survival)
## (gSurv.B2 <- clogit(y ~ x + strata(block), data=dB2))
## ## --> Error in Surv(rep(1, 200L), y) : Time and status are different lengths
## summary(gSurv.B2)
## (SE.surf <- sqrt(diag(vcov(gSurv.B2))))



g3 <-  glmer(y ~ x + (1|block),data=gSim(family=binomial(),nbinom=10),
             family=binomial)
expect_equal(var(sumFun(g3)),0)

## check dispersion parameter
## (lowered tolerance to pass checks on my machine -- SCW)
expect_equal(sigma(g0)^2, 0.4888248, tolerance=1e-4)

