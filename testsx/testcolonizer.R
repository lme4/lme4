library(lme4.0)
## Emacs M-<Enter> --> setwd() correctly
(load("colonizer_rand.RData"))
m0.0 <- glm(colonizers~Treatment*homespecies*respspecies, data=randdat, family=poisson)
with(randdat,tapply(colonizers,list(Treatment,homespecies,respspecies),sum))
summary(m1.0 <- glmer(form1, data=randdat, family=poisson))
summary(m2.0 <- glmer(form2, data=randdat, family=poisson))

detach("package:lme4.0", unload=TRUE)
library(lme4)
packageDescription("lme4")

## currently (r 1704), all give  "Downdated VtV is not positive definite"
try(m1 <- glmer(form1,data=randdat, family=poisson, verbose=10L))
try(m1 <- glmer(form1,data=randdat, family=poisson, verbose=10L, tolPwrss=1e-13))
try(m2 <- glmer(form2,data=randdat, family=poisson, verbose=10L))
detach("package:lme4", unload=TRUE)
