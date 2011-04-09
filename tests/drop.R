library(lme4)
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)

gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
             family = binomial, data = cbpp)

terms(fm1)
extractAIC(fm1)
drop1(fm1)
drop1(fm1,test="Chisq")

drop1(gm1,test="Chisq")
