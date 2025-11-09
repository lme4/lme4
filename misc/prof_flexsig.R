devtools::load_all()

data("sleepstudy", package = "lme4")
m1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
reCovs <- getReCovs(m1)
getProfPars(m1, profscale = "sdcor")
p <- getProfPars(reCovs[[1]], profscale = "sdcor")
setProfPars(reCovs[[1]], p, profscale = "sdcor")
setProfPars(m1, c(p, sigma(m1)), profscale = "sdcor")

d2 <- devfun2(m1)
## np == 4
np <- environment(d2)$np
p0 <- unlist(attr(d2, "optimum")[1:np])
trace("setProfPars")
d2(p0)
f <- function(x) d2(p0 - c(x, 0, 0,0))
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


