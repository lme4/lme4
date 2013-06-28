## see if we can still run lme4 functions when lme4 is not attached
if ("package:lme4" %in% search()) detach("package:lme4")
data(sleepstudy,package="lme4")
data(cbpp,package="lme4")
fm1 <- lme4::lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
gm1 <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
             data = cbpp, family = binomial)
