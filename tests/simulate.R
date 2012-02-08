library(lme4Eigen)
## debug(simulate.merMod)
gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
              data = cbpp, family = binomial)
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)

simulate(gm1)
simulate(fm1)

## FIXME: would like better tests ...
