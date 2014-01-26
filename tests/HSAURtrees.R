library("lme4")

## example from HSAUR2 package; data from 'multcomp'; see ../inst/testdata/trees513.R
load(system.file("testdata","trees513.RData",package="lme4"))

## model formula:
modForm <- damage ~ species - 1 + (1 | lattice / plot)

dfun <- glmer(modForm, data = trees513B, family = binomial,
              devFunOnly = TRUE)
ls.str(environment(dfun))# "for your information"

.not.call <- function(x) x[names(x) != "call"]

if(lme4:::testLevel() < 2) q("no")
## {{advantage to  if(. >= 2) { ........} : autoprint of system.time() etc

## else  (testLevel >= 2) : --------------------------------------------------

## Generate oldres:
## ----------------
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
system.time(mmodA <- glmer(modForm, data = trees513A, family = binomial()))
## 7 seconds
newres <- c(fixef(mmodA), getME(mmodA,"theta"))
stopifnot(all.equal(oldres, newres, tolerance=1.5e-3))
system.time(mmodB <- glmer(modForm, data = trees513B, family = binomial()))
## 10.4 seconds
##
## lmer( + family) -> diverts to glmer() with a warning [TODO: use assertWarning(.) eventually]
system.time(lmodB <-
            lmer(modForm, data = trees513B, family = binomial()))
stopifnot(all.equal(.not.call(summary(mmodB)),
                    .not.call(summary(lmodB))))
newresB <- c(fixef(mmodB),getME(mmodB,"theta"))
stopifnot(length(newresB) == length(oldres) + 1)# extra: species[ash/maple/elm/lime]
## (unfinished)


