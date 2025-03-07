library(lme4)
## show important current settings {for reference, etc} -- [early, and also on Windows !]:

source(system.file("test-tools-1.R", package = "Matrix"), keep.source = FALSE)
## N.B. is.all.equal[34]() and assert.EQ() use 'tol', not 'tolerance'

str( lmerControl())
str(glmerControl())
str(nlmerControl())
ls.str(environment(nloptwrap))


## see Details under ?deviance.merMod:

fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
fm1ML <- refitML(fm1)
REMLcrit(fm1)   ## 1743.628
deviance(fm1ML) ## 1751.939
deviance(fm1,REML=FALSE)  ## FIXME: not working yet (NA)
deviance(fm1,REML=TRUE)

## from lme4.0
oldvals <- c(REML=1743.6282722424, ML=1751.98581103058)
## leave out ML values for REML fits for now ...
stopifnot(
          is.all.equal3(REMLcrit(fm1), deviance(fm1,REML=TRUE), deviance(fm1ML,REML=TRUE),oldvals["REML"]),
          all.equal(deviance(fm1ML),deviance(fm1ML,REML=FALSE),oldvals["ML"]),
          all.equal(REMLcrit(fm1)/-2,c(logLik(fm1)),c(logLik(fm1ML,REML=TRUE)),c(logLik(fm1,REML=TRUE))),
          all.equal(deviance(fm1ML)/-2,c(logLik(fm1ML,REML=FALSE)),
                    c(logLik(fm1ML,REML=FALSE))))

## should be:
## stopifnot(
##           all.equal(deviance(fm1),deviance(fm1,REML=TRUE),deviance(fm1ML,REML=TRUE),oldvals["REML"]),
##           all.equal(deviance(fm1ML),deviance(fm1,REML=FALSE),deviance(fm1ML,REML=FALSE),oldvals["ML"]),
##           all.equal(deviance(fm1)/2,c(logLik(fm1)),c(logLik(fm1ML,REML=TRUE)),c(logLik(fm1,REML=TRUE))),
##           all.equal(deviance(fm1ML)/2,c(logLik(fm1,REML=FALSE)),c(logLik(fm1ML,REML=FALSE)),
##                     c(logLik(fm1ML,REML=FALSE))))
