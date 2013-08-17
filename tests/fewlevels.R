## example from Gabor Grothendieck; didn't 
## https://stat.ethz.ch/pipermail/r-sig-mixed-models/2010q2/003726.html
## library(nlme)
## set.seed(1)
## f <- function(n, k) {
##    set.seed(1)
##    x <- 1:n
##    fac <- gl(k, 1, n)
##    fac.eff <- rnorm(k, 0, 4)[fac]
##    e <- rnorm(n)
##    y <- 1 + 2 * x + fac.eff + e
##    lme(y ~ x, random = ~ 1 | fac)
## }
## n <- 10000
sim <- function(n=1e4,k=4, seed=1) {
    set.seed(seed)
    x <- 1:n
    fac <- gl(k, 1, n)
    fac.eff <- rnorm(k, 0, 4)[fac]
    e <- rnorm(n)
    y <- 1 + 2 * x + fac.eff + e
    data.frame(x,y,fac)
}
library(nlme)
m.lme <- lme(y ~ x, random=~ 1|fac ,data=sim())
v.lme <- as.numeric(VarCorr(m.lme)[1,1])
library(lme4)

m.lmer <- lmer(y ~ x + (1|fac),data=sim())
v.lmer <- VarCorr(m.lmer)[[1]][1,1]
stopifnot(all.equal(v.lmer,19.54829,tol=1e-6))
