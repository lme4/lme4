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

m1 <- glmer(form1,data=randdat, family=poisson)  ## PIRLS step failed
m2 <- glmer(form1,data=randdat, family=poisson, nAGQ=0)  ## OK
m3 <- glmer(form2,data=randdat, family=poisson)  ## ditto
