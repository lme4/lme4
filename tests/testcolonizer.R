## library(lme4.0)
## Emacs M-<Enter> --> setwd() correctly

## m0.0 <- glm(colonizers~Treatment*homespecies*respspecies, data=randdat, family=poisson)
## with(randdat,tapply(colonizers,list(Treatment,homespecies,respspecies),sum))
## summary(m1.0 <- glmer(form1, data=randdat, family=poisson))
## summary(m2.0 <- glmer(form2, data=randdat, family=poisson))


## detach("package:lme4.0", unload=TRUE)

load(system.file("testdata","colonizer_rand.rda",package="lme4"))
library("lme4")
packageVersion("lme4")

## FIXME: currently fails.  Works in lme4.0, BUT ... this is a very poorly
##   posed problem (complete separation etc.)  nAGQ=0 gives decent results.
try(m1 <- glmer(form1,data=randdat, family=poisson, verbose=10L))  ## PIRLS step failed
try(m1 <- glmer(form1,data=randdat, family=poisson, verbose=10L, nAGQ=0))  ## OK
try(m2 <- glmer(form2,data=randdat, family=poisson, verbose=10L))  ## ditto
