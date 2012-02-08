library("testthat")
context("glmFamily")
eps <- .Machine$double.eps
oneMeps <- 1 - eps
set.seed(1)

## sample linear predictor values for the unconstrained families
etas <- list(seq.int(-8, 8, by=1),  # equal spacing to asymptotic area
             runif(17, -8, 8),  # random sample from wide uniform dist
             rnorm(17, 0, 8),   # random sample from wide normal dist
             c(-10^30, rnorm(15, 0, 4), 10^30))

## sample linear predictor values for the families in which eta must be positive
etapos <- list(seq.int(1, 20, by=1),
               rexp(20),
               rgamma(20, 3),
               pmax(.Machine$double.eps, rnorm(20, 2, 1)))

## values of mu in the (0,1) interval
mubinom <- list(runif(100, 0, 1),
                rbeta(100, 1, 3),
                pmin(pmax(eps, rbeta(100, 0.1, 3)), oneMeps),
                pmin(pmax(eps, rbeta(100, 3, 0.1)), oneMeps))

test_that("inverse link and muEta functions", {
    tst.lnki <- function(fam, frm) {
        ff <- glmFamily$new(family=fam)
        sapply(frm, function(x) expect_that(fam$linkinv(x), equals(ff$linkInv(x))))
    }

    tst.muEta <- function(fam, frm) {
        ff <- glmFamily$new(family=fam)
        sapply(frm, function(x) expect_that(fam$mu.eta(x), equals(ff$muEta(x))))
    }
    
    tst.lnki(binomial(),           etas) # binomial with logit link
    tst.muEta(binomial(),          etas)
    tst.lnki(binomial("probit"),   etas) # binomial with probit link
    tst.muEta(binomial("probit"),  etas)
    tst.lnki(binomial("cloglog"),  etas) # binomial with cloglog link
    tst.muEta(binomial("cloglog"), etas)
    tst.lnki(binomial("cauchit"),  etas) # binomial with cauchit link
    tst.muEta(binomial("cauchit"), etas)
    tst.lnki(poisson(),            etas) # Poisson with log link
    tst.muEta(poisson(),           etas)
    tst.lnki(gaussian(),           etas) # Gaussian with identity link
    tst.muEta(gaussian(),          etas)
    tst.lnki(Gamma(),              etapos) # gamma family
    tst.muEta(Gamma(),             etapos)
    tst.lnki(inverse.gaussian(),   etapos) # inverse Gaussian
    tst.muEta(inverse.gaussian(),  etapos)    
})

test_that("link and variance functions", {
    tst.link <- function(fam, frm) {
        ff <- glmFamily$new(family=fam)
        sapply(frm, function(x) expect_that(fam$linkfun(x), equals(ff$link(x))))
    }

    tst.variance <- function(fam, frm) {
        ff <- glmFamily$new(family=fam)
        sapply(frm, function(x) expect_that(fam$variance(x), equals(ff$variance(x))))
    }

    tst.link(    binomial(),          mubinom)
    tst.variance(binomial(),          mubinom)
    tst.link(    binomial("probit"),  mubinom)
    tst.variance(binomial("probit"),  mubinom)
    tst.link(    binomial("cauchit"), mubinom)
    tst.variance(binomial("cauchit"), mubinom)
    tst.link(    gaussian(),          etas)
    tst.variance(gaussian(),          etas)
    tst.link(    Gamma(),             etapos)
    tst.variance(Gamma(),             etapos)
    tst.link(    inverse.gaussian(),  etapos)
    tst.variance(inverse.gaussian(),  etapos)    
})



