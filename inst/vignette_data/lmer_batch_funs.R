library("lme4")

fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

ccw <- confint(fm1, method = "Wald")
## compute the profile once here and re-use it below rather than
## profiling fm1 a second time; profiling is expensive and
## makes us exceed CRAN time limits
pf <- profile(fm1)
ccp <- confint(pf)
ccb <- confint(fm1, method = "boot")

save(ccw, pf, ccp, ccb, file = "lmer_batch.rda")

sessionInfo()
