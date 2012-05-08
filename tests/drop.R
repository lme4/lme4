library(lme4)
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)

## slightly weird model but plausible --- not that
##   one would want to try drop1() on this model ...
fm2 <- lmer(Reaction ~ 1+ (Days|Subject), sleepstudy)
drop1(fm2)  ## empty
update(fm1, . ~ . - Days)
anova(fm2) ## empty

terms(fm1)
terms(fm1,fixed.only=FALSE)

extractAIC(fm1)

drop1(fm1)
drop1(fm1, test="Chisq")

gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
             family = binomial, data = cbpp, nAGQ=25L)

drop1(gm1, test="Chisq")

