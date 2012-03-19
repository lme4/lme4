library(lme4)
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
stopifnot(isREML(fm1),
          isLMM(fm1),
          !isGLMM(fm1),
          !isNLMM(fm1))

fm1ML <- refitML(fm1)
stopifnot(!isREML(fm1ML),
          isLMM(fm1ML),
          !isGLMM(fm1ML),
          !isNLMM(fm1ML))

gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
              data = cbpp, family = binomial)
stopifnot(!isREML(gm1),
          !isLMM(gm1),
          isGLMM(gm1),
          !isNLMM(gm1))

nm1 <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree,
             Orange, start = c(Asym = 200, xmid = 725, scal = 350))
stopifnot(!isREML(nm1),
          !isLMM(nm1),
          !isGLMM(nm1),
          isNLMM(nm1))


