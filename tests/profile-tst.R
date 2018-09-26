library(lme4)
library(testthat)
library(lattice)

options(nwarnings = 5000)# instead of 50, and then use  summary(warnings())

### __ was ./profile_plots.R ___
fm1 <- lmer(Reaction~ Days + (Days|Subject), sleepstudy)
pfile <- system.file("testdata","tprfm1.RData", package="lme4")
if(file.exists(pfile)) print(load(pfile)) else withAutoprint({
    system.time( tpr.fm1 <- profile(fm1, optimizer="Nelder_Mead") ) ## 5 sec (2018); >= 50 warnings !?
    save(tpr.fm1, file= "../../inst/testdata/tprfm1.RData")
})
oo <- options(warn = 2) # {warnings are errors from here on}

if(!dev.interactive(orNone=TRUE)) pdf("profile_plots.pdf")
xyplot(tpr.fm1)
splom(tpr.fm1)
densityplot(tpr.fm1, main="densityplot( profile(lmer(..)) )")

## various scale options
xyplot(tpr.fm1,scale=list(x=list(relation="same")))  ## stupid
xyplot(tpr.fm1,scale=list(y=list(relation="same")))
xyplot(tpr.fm1,scale=list(y=list(relation="same"),tck=0))

##
expect_error(xyplot(tpr.fm1,conf=50),"must be strictly between 0 and 1")

### end {profile_plots.R}

fm01ML <- lmer(Yield ~ 1|Batch, Dyestuff, REML = FALSE)

## 0.8s (on a 5600 MIPS 64bit fast(year 2009) desktop "AMD Phenom(tm) II X4 925"):
##
system.time( tpr <- profile(fm01ML) )

## test all combinations of 'which', including plots (but don't show plots)
wlist <- list(1:3,1:2,1,2:3,2,3,c(1,3))
invisible(lapply(wlist,function(w) xyplot(profile(fm01ML,which=w))))

(confint(tpr) -> CIpr)
print(xyplot(tpr))
##  comparing against lme4a reference values -- but lme4 returns sigma
## rather than log(sigma)
stopifnot(dim(CIpr) == c(3,2),
          all.equal(unname(CIpr[".sigma",]),exp(c(3.64362, 4.21446)), tolerance=1e-6),
          all.equal(unname(CIpr["(Intercept)",]),c(1486.451500,1568.548494)))

options(oo)# warnings allowed ..

## fixed-effect profiling with vector RE
data(Pastes)
fmoB <- lmer(strength ~ 1 + (cask | batch), data=Pastes,
             control = lmerControl(optimizer = "bobyqa"))
(pfmoB <- profile(fmoB, which = "beta_", alphamax=.001))
xyplot(pfmoB)# nice and easy ..

summary(
    fm <- lmer(strength ~ 1 + (cask | batch), data=Pastes,
               control = lmerControl(optimizer = "nloptwrap",
                                     calc.derivs= FALSE))
)

ls.str(environment(nloptwrap))# showing *its* defaults

pfm <- profile(fm, which = "beta_", alphamax=.001) # 197 warnings for "nloptwrap"
summary(warnings())
str(pfm) # only 3 rows, .zeta = c(0, NaN, Inf) !!!
try( xyplot(pfm) ) ## FIXME or rather the profiling or rather the "wrap on nloptr"

(testLevel <- lme4:::testLevel())
if(testLevel > 2) {

    ## 2D profiles
    fm2ML <- lmer(diameter ~ 1 + (1|plate) + (1|sample), Penicillin, REML=0)
    system.time(pr2 <- profile(fm2ML)) # 5.2 sec, 2018-05: 2.1"
    (confint(pr2) -> CIpr2)

    lme4a_CIpr2 <-
        structure(c(0.633565787613112, 1.09578224011285, -0.721864513060904,
                    21.2666273835452, 1.1821039843372, 3.55631937954106, -0.462903300019305,
                    24.6778176174587), .Dim = c(4L, 2L), .Dimnames = list(c(".sig01",
                           ".sig02", ".lsig", "(Intercept)"), c("2.5 %", "97.5 %")))
    lme4a_CIpr2[".lsig",] <- exp(lme4a_CIpr2[".lsig",])

    stopifnot(all.equal(unname(CIpr2),unname(lme4a_CIpr2),tolerance=1e-6))

    print(xyplot(pr2, absVal=0, aspect=1.3, layout=c(4,1)))
    print(splom(pr2))

    gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                 data = cbpp, family = binomial)

    ## GLMM profiles
    system.time(pr4 <- profile(gm1))  ## ~ 10 seconds

    pr4.3 <- profile(gm1,which=3)
    xyplot(pr4,layout=c(5,1),as.table=TRUE)

    splom(pr4) ## used to fail because of NAs

    nm1 <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree,
                 Orange, start = c(Asym = 200, xmid = 725, scal = 350))
    if (FALSE) {
        ## not working yet: detecting (slightly) lower deviance; not converging in 10k
        pr5 <- profile(nm1,which=1,verbose=1,maxmult=1.2)
        xyplot(.zeta~.focal|.par,type=c("l","p"),data=lme4:::as.data.frame.thpr(pr5),
               scale=list(x=list(relation="free")),
               as.table=TRUE)
    }
}  ## testLevel > 2

if (testLevel > 3) {
    fm3ML <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, REML=FALSE)
    ## ~ 4 theta-variables (+ 2 fixed), 19 seconds | 2018-05: 7.4"
    print(system.time(pr3 <- profile(fm3ML)))
    print(xyplot(pr3))
    print(splom(pr3))

    if (testLevel > 4) {
      if(requireNamespace("mlmRev")) {
        data("Contraception", package="mlmRev")
        ## fit already takes ~ 3 sec (2018-05)
        fm2 <- glmer(use ~ urban+age+livch + (urban|district), Contraception, binomial)
        print(system.time(pr5 <- profile(fm2,verbose=10))) # 2018-05: 462 sec = 7'42"
        ## -> 5 warnings notably "non-monotonic profile for .sig02" (the RE's corr.)
        print(xyplot(pr5))
      }
    }  ## testLevel > 4

}  ## testLevel > 3

library("parallel")
if (detectCores()>1) {

    p0 <- profile(fm1, which="theta_")
    ## http://stackoverflow.com/questions/12983137/how-do-detect-if-travis-ci-or-not
    travis <- nchar(Sys.getenv("TRAVIS")) > 0
    if(.Platform$OS.type != "windows" && !travis) {
        prof01P <- profile(fm1, which="theta_", parallel="multicore", ncpus=2)
        stopifnot(all.equal(p0,prof01P))
    }

    ## works in Solaris from an interactive console but not ???
    ##   via R CMD BATCH

    if (Sys.info()["sysname"] != "SunOS" && !travis) {
        prof01P.snow <- profile(fm1, which="theta_", parallel="snow", ncpus=2)
        stopifnot(all.equal(p0,prof01P.snow))
    }
}

## test profile/update from within functions
foo <- function() {
    gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                 data = cbpp, family = binomial)
    ## return
    profile(gm1, which="theta_")
}
stopifnot(inherits(foo(), "thpr"))
