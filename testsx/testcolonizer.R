library(lme4.0)
load("colonizer_rand.RData")
m1 <- glmer(form1,data=randdat,
      family=poisson)
m2 <- glmer(form2,data=randdat,
      family=poisson)
detach("package:lme4.0")
library(lme4)
try(m1 <- glmer(form1,data=randdat,
      family=poisson,verbose=10L))
try(m1 <- glmer(form1,data=randdat,
      family=poisson,verbose=10L,tolPwrss=1e-13))
try(m2 <- glmer(form2,data=randdat,
      family=poisson,verbose=10L))

