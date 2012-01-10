library(lme4Eigen)

fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)

terms(fm1)
extractAIC(fm1)
drop1(fm1)
drop1(fm1, test="Chisq")


gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
             family = binomial, data = cbpp, nAGQ=25L)

drop1(gm1, test="Chisq")
