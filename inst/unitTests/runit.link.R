eps <- .Machine$double.eps
oneMeps <- 1 - eps
set.seed(1)
etas <-
    lapply(list(-8:8,             # equal spacing to asymptotic area
                runif(20, -8, 8), # random sample from wide uniform dist
                rnorm(20, 0, 8),  # random sample from wide normal dist
                c(-10^30, rnorm(10, 0, 4), 10^30)), as.numeric)
etapos <-
    lapply(list(1:20,
                rexp(20),
                rgamma(20, 3),
                pmax(eps, rnorm(20, 2, 1))), as.numeric)

mubinom <-
    lapply(list(runif(100, 0, 1),
                rbeta(100, 1, 3),
                pmin(pmax(eps, rbeta(100, 0.1, 3)), oneMeps),
                pmin(pmax(eps, rbeta(100, 3, 0.1)), oneMeps)), as.numeric)

tst.lnki <- function(fam, lst) {
    ff <- glmFamily$new(fam)
    unlist(lapply(lst, function(x)
                  checkEquals(fam$linkinv(x), ff$linkInv(x))))
}
tst.link <- function(fam, lst) {
    ff <- glmFamily$new(fam)
    unlist(lapply(lst, function(x)
                  checkEquals(fam$linkfun(x), ff$linkFun(x))))
}
tst.muEta <- function(fam, lst) {
    ff <- glmFamily$new(fam)
    unlist(lapply(lst, function(x)
                  checkEquals(fam$mu.eta(x), ff$muEta(x))))
}

test.uncons.lnki <- function() {        # linkinv on unconstrained eta
    tst.lnki(binomial(), etas)          # binomial with logit link
    tst.muEta(binomial(), etas)
    tst.lnki(binomial("probit"), etas)  # binomial with probit link
    tst.muEta(binomial("probit"), etas)
    tst.lnki(binomial("cloglog"), etas) # binomial with cloglog link
    tst.muEta(binomial("cloglog"), etas)
    tst.lnki(binomial("cauchit"), etas) # binomial with cauchit link
    tst.muEta(binomial("cauchit"), etas)
    tst.lnki(poisson(), etas)           # Poisson with log link
    tst.muEta(poisson(), etas)
    tst.lnki(gaussian(), etas)          # Gaussian with identity link
    tst.muEta(gaussian(), etas)
}

test.pos.lnki <- function() {           # linkinv for positive eta
    set.seed(1)
    tst.lnki(Gamma(), etapos)           # gamma family
    tst.muEta(Gamma(), etapos)
    tst.lnki(inverse.gaussian(), etapos) # inverse Gaussian
    tst.muEta(inverse.gaussian(), etapos)    
}

test.binom.link <- function() {         # link for binomial mu
    tst.link(binomial(), mubinom)
    tst.link(binomial("probit"), mubinom)
}

test.pos.link <- function() {           # link for positive mu
    tst.link(poisson(), etapos)
    tst.link(Gamma(), etapos)
    tst.link(inverse.gaussian(), etapos)    
}

test.uncons.link <- function() {        # link for unconstrained mu
    tst.link(gaussian(), etas)
}

## ToDo: Add checks on variance functions
