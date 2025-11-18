devtools::load_all()

data("sleepstudy", package = "lme4")
m0 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
p0 <- profile(m0)
m1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
p1 <- profile(m1)
m2 <- lmer(Reaction ~ Days + (Days || Subject), sleepstudy)
p2 <- profile(m2)
m3 <- lmer(Reaction ~ Days + cs(Days | Subject), sleepstudy)
p3 <- profile(m3)


dd <- devfun2(m1)
par1 <- c(16.066, 4.744, 0.955, 27.525)
mkpar <- function(rho) { par1[3] <- rho; par1 }
dd(par1)
dd(mkpar(0.999))
dd(mkpar(1.0)) ## NA

xyplot(p1)
debug(profile)

p2 <- profile(m2)

master_ci <- structure(c(14.3814181607732, -0.481500758280773, 3.8011640532591, 
22.8982668820808, 237.680695505234, 7.35865327511699, 37.7159953182866, 
0.684986273295485, 8.75338075305527, 28.8579965114253, 265.129514680174, 
13.5759187517082), dim = c(6L, 2L), dimnames = list(c(".sig01", 
".sig02", ".sig03", ".sigma", "(Intercept)", "Days"), c("2.5 %", 
                                                        "97.5 %")))

mm <- master_ci[c(1,3,2,4:6),]
abs((mm-ci1)/mm)  ## relative differences (why?)

ci1 <- confint(m1)
all.equal(ci1, ,
          check.attributes = FALSE,
          tolerance = 2e-6)

profile(m0, ".sigma")
lattice::xyplot(p0)
reCovs <- getReCovs(m2)
getProfPars(m1, profscale = "sdcor")
getProfPars(m2, profscale = "sdcor")
getProfPars(reCovs[[1]], profscale = "sdcor")
p <- getProfPars(reCovs[[1]], profscale = "sdcor")
setProfPars(reCovs[[1]], p, profscale = "sdcor")
setProfPars(m1, c(p, sigma(m1)), profscale = "sdcor")

d2 <- devfun2(m1)
d3 <- devfun2(m2)
## np == 4
np <- environment(d3)$np
p0 <- unlist(attr(d3, "optimum")[1:np])
trace("setProfPars")
d3(p0)

np <- environment(d3)$np
p0 <- unlist(attr(d2, "optimum")[1:np])
trace("setProfPars")
d2(p0)

f <- function(x) d3(p0 - c(x, 0, 0,0))
xvec <- seq(-30, 20, length = 51)
lvec <- sapply(xvec, f)
plot(xvec, lvec)

trace("setVC", sig = c("Covariance.us", "numeric", "numeric"), browser)
untrace("setVC", sig = c("Covariance.us", "numeric", "numeric"))

debug(profile.merMod)
p1 <- profile(m1, which = 1, verbose = TRUE)
## lots of warnings but ?? maybe there all along ???
lattice::xyplot(p1)
confint(p1)
##         2.5 %   97.5 %
## .sig01 14.3818 37.71583
p2 <- profile(m1, which = 2, verbose = TRUE)
lattice::xyplot(p2)
p3 <- profile(m1, which = 3, verbose = TRUE)
lattice::xyplot(p3)

p_all <- profile(m1, verbose = TRUE)
confint(p_all)
xyplot(p_all)

v <- c(23.7805585271338, 5.71683789263024, 0.0813195978647574)
s <- 25.59182
m <- diag(2)
m[2,1] <- m[1,2] <- v[3]
vv <- v[1:2]/s
chol(m * outer(vv, vv))
c(0.929, 0.018, 0.222)
profile(m1)



library(testthat)
library(lme4)
fm1 <- lmer(Reaction ~ Days + (1|Subject), sleepstudy)
p1 <- lme4:::getPar(fm1)
p2 <- lme4:::convParToProfPar(p1, fm1, "sdcor", sigma(fm1))
p3 <- lme4:::convProfParToPar(p2, fm1, "sdcor", sigma(fm1))
stopifnot(all.equal(p1,p3))

devtools::load_all()
debug(profile.merMod)
profile(fm1)
