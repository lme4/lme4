library(lme4)
library(lattice)

### __ was ./profile_plots.R ___
fm1 <- lmer(Reaction~ Days + (Days|Subject), sleepstudy)
pfile <- system.file("testdata","tprfm1.RData", package="lme4")
if(file.exists(pfile)) print(load(pfile)) else {
 system.time( tpr.fm1 <- profile(fm1, optimizer="Nelder_Mead") )  ## 20 seconds
 save(tpr.fm1, file= "../../inst/testdata/tprfm1.RData")
}
oo <- options(warn = 1) # {warnings are errors from here on}
                        # FIXME: switched warnings back to get through checks

if(!dev.interactive(orNone=TRUE)) pdf("profile_plots.pdf")
xyplot(tpr.fm1)
splom(tpr.fm1)
densityplot(tpr.fm1, main="densityplot( profile(lmer(..)) )")
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
fm <- lmer(strength ~ 1 + (cask | batch), data=Pastes)
pfm <- profile(fm, which = "beta_", alphamax=.001)
xyplot(pfm)

(testLevel <- lme4:::testLevel())
if(testLevel > 2) {

    ## 2D profiles
    fm2ML <- lmer(diameter ~ 1 + (1|plate) + (1|sample), Penicillin, REML=0)
    system.time(pr2 <- profile(fm2ML)) # 5.2 sec
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
        ## not working yet
        pr5 <- profile(nm1,which=1,verbose=1,maxmult=1.2)
        xyplot(.zeta~.focal|.par,type=c("l","p"),data=lme4:::as.data.frame.thpr(pr5),
               scale=list(x=list(relation="free")),
               as.table=TRUE)
    }
}  ## testLevel > 2

if (testLevel > 3) {
    fm3ML <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, REML=FALSE)
    ## ~ 4 theta-variables (+ 2 fixed), 19 seconds :
    print(system.time(pr3 <- profile(fm3ML)))
    print(xyplot(pr3))
    print(splom(pr3))

    if (testLevel > 4) {
        if(require("mlmRev")) {
            ## takes much longer
            data("Contraception", package="mlmRev")
            fm2 <- glmer(use ~ urban+age+livch+(urban|district), Contraception, binomial)
            print(system.time(pr5 <- profile(fm2,verbose=10)))
            print(xyplot(pr5))
        }
    }  ## testLevel > 4

}  ## testLevel > 3
