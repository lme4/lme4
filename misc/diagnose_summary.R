Sys.setenv("LME4_TEST_LEVEL" = 100)
devtools::test()

library(lme4)
m1 <- lmer(Reaction ~ Days+(Days|Subject), sleepstudy)
debug(lme4:::summary.merMod)
summary(m1)
cc <- capture.output(print(summary(m1)))
any(grepl("Correlation of Fixed Effects", cc))

library(merDeriv)
cc <- capture.output(print(summary(m1)))
any(grepl("Correlation of Fixed Effects", cc))

debug(lme4:::print.summary.merMod)
debug(lme4:::vcov.merMod)
length(vcov(m1)@factors)
s1 <- summary(m1)
str(s1$vcov)  ## no factors when run AFTER test() at testLevel 100
           
f <- function(L) {
    load(sprintf("test-summary_testlevel_%d.rda", L))
    return(mget(ls()))
}

r1 <- f(1)
r100 <- f(100)
waldo::compare(r1$m1, r100$m1)
summary(r1$m1)
summary(r100$m1)

## cc1 doesn't include 

waldo::compare(r1$ss, r100$ss)
waldo::compare(r1$m1, r100$m1, ignore_attr = TRUE)

setdiff(names(r100$ss$loadedOnly), names(r1$ss$loadedOnly))  ## MASS, boot, splines
setdiff(names(r100$ss$otherPkgs), names(r1$ss$otherPkgs))  ## MASS, boot, splines
## what is messing things up? stats4, gamm4?
