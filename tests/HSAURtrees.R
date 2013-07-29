library("lme4")
(testLevel <- if (nzchar(s <- Sys.getenv("LME4_TEST_LEVEL"))) as.numeric(s) else 1)

## example from HSAUR2 package; data from 'multcomp'; see ../inst/testdata/trees513.R
load(system.file("testdata","trees513.RData",package="lme4"))

## dfun <- glmer(damage ~ species - 1 + (1 | lattice / plot),
##              data = trees513, family = binomial, devFunOnly = TRUE)
## ls.str(environment(dfun))# and you can investigate ...

if (testLevel > 2) {
    ## library(lme4.0)
    ## system.time(mmod0 <- glmer(damage ~ species - 1 + (1 | lattice / plot),
    ##               data = trees513, family = binomial()))
    ## ## 4 seconds
    ## oldres <- c(fixef(mmod0),getME(mmod0,"theta"))
    ## detach("package:lme4.0")
    ## dput(oldres)
    oldres <- structure(c(5.23645064474105, 4.73568475545248, 2.65289926317093,
                          1.29043984816924, 1.59329381563025,
                          0.532663142106669, 1.16703186884403
                          ), .Names = c("speciesspruce", "speciespine",
                             "speciesbeech",
                             "speciesoak", "specieshardwood",
                             "plot:lattice.(Intercept)",
                             "lattice.(Intercept)"))
    system.time(mmodA <- glmer(damage ~ species - 1 + (1 | lattice / plot),
                              data = trees513A, family = binomial()))
    ## 7 seconds
    newres <- c(fixef(mmodA),getME(mmodA,"theta"))
    stopifnot(all.equal(oldres,newres,tol=1.5e-3))
    system.time(mmodB <- glmer(damage ~ species - 1 + (1 | lattice / plot),
                              data = trees513B, family = binomial()))
    ## 10.4 seconds
    newresB <- c(fixef(mmodB),getME(mmodB,"theta"))
}

