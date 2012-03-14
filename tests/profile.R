library(lme4Eigen)

fm01ML <- lmer(Yield ~ 1|Batch, Dyestuff, REML = FALSE)

## 0.8s (on a 5600 MIPS 64bit fast(year 2009) desktop "AMD Phenom(tm) II X4 925"):
## 
system.time( tpr <- profile(fm01ML) )

(confint(tpr) -> CIpr)
print(xyplot(tpr))
##  comparing against lme4a reference values -- but lme4Eigen returns sigma
## rather than log(sigma)
stopifnot(dim(CIpr) == c(3,2),
          all.equal(unname(CIpr[".sigma",]),exp(c(3.64362, 4.21446)), tol=1e-6),
          all.equal(unname(CIpr["(Intercept)",]),c(1486.451500,1568.548494)))

## 2D profiles
fm2ML <- lmer(diameter ~ 1 + (1|plate) + (1|sample), Penicillin, REML=0)
pr2 <- profile(fm2ML)
(confint(pr2) -> CIpr2)

lme4a_CIpr2 <-
structure(c(0.633565787613112, 1.09578224011285, -0.721864513060904, 
21.2666273835452, 1.1821039843372, 3.55631937954106, -0.462903300019305, 
24.6778176174587), .Dim = c(4L, 2L), .Dimnames = list(c(".sig01", 
".sig02", ".lsig", "(Intercept)"), c("2.5 %", "97.5 %")))
lme4a_CIpr2[".lsig",] <- exp(lme4a_CIpr2[".lsig",])

stopifnot(all.equal(unname(CIpr2),unname(lme4a_CIpr2)))

print(xyplot(pr2, absVal=0, aspect=1.3, layout=c(4,1)))
print(splom(pr2))

## NOT RUN: takes ~ 30 seconds on my machine ...
fm3ML <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, REML=FALSE)
## system.time(pr3 <- profile(fm3ML))
## xyplot(pr3)
## print(splom(pr3))

