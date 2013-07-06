library(lme4)
(testLevel <- if (nzchar(s <- Sys.getenv("LME4_TEST_LEVEL"))) as.numeric(s) else 1)

fm01ML <- lmer(Yield ~ 1|Batch, Dyestuff, REML = FALSE)

## 0.8s (on a 5600 MIPS 64bit fast(year 2009) desktop "AMD Phenom(tm) II X4 925"):
##
system.time( tpr <- profile(fm01ML) )

## test all combinations of 'which', including plots
wlist <- list(1:3,1:2,1,2:3,2,3,c(1,3))
invisible(lapply(wlist,function(w) xyplot(profile(fm01ML,which=w))))

(confint(tpr) -> CIpr)
print(xyplot(tpr))
##  comparing against lme4a reference values -- but lme4 returns sigma
## rather than log(sigma)
stopifnot(dim(CIpr) == c(3,2),
          all.equal(unname(CIpr[".sigma",]),exp(c(3.64362, 4.21446)), tol=1e-6),
          all.equal(unname(CIpr["(Intercept)",]),c(1486.451500,1568.548494)))

if (testLevel > 2) {

    ## 2D profiles
    fm2ML <- lmer(diameter ~ 1 + (1|plate) + (1|sample), Penicillin, REML=0)
    system.time(pr2 <- profile(fm2ML))
    (confint(pr2) -> CIpr2)

    lme4a_CIpr2 <-
        structure(c(0.633565787613112, 1.09578224011285, -0.721864513060904,
                    21.2666273835452, 1.1821039843372, 3.55631937954106, -0.462903300019305,
                    24.6778176174587), .Dim = c(4L, 2L), .Dimnames = list(c(".sig01",
                           ".sig02", ".lsig", "(Intercept)"), c("2.5 %", "97.5 %")))
    lme4a_CIpr2[".lsig",] <- exp(lme4a_CIpr2[".lsig",])

    stopifnot(all.equal(unname(CIpr2),unname(lme4a_CIpr2),tol=1e-6))

    print(xyplot(pr2, absVal=0, aspect=1.3, layout=c(4,1)))
    print(splom(pr2))

    gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                 data = cbpp, family = binomial)

    ## GLMM profiles
    system.time(pr4 <- profile(gm1))  ## ~ 10 seconds

    ## FIXME: compDev=FALSE fails
    if (FALSE) {
        gm1B <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                      data = cbpp, family = binomial,
                      compDev=FALSE,verbose=3)
    }

    profile(gm1,which=3)
    xyplot(pr4,layout=c(5,1),as.table=TRUE)
    if (FALSE) {
        ## FIXME! fails because of NAs
        splom(pr4)
    }

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

## NOT RUN:  ~ 4 theta-variables, 19 seconds
fm3ML <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, REML=FALSE)
if (testLevel > 3) {
    system.time(pr3 <- profile(fm3ML))
    xyplot(pr3)
    print(splom(pr3))

    data("Contraception",package="mlmRev")
    fm2 <- glmer(use ~ urban+age+livch+(urban|district), Contraception, binomial)
    pr5 <- profile(fm2,verbose=10)
    xyplot(pr5)
}  ## testLevel > 3
