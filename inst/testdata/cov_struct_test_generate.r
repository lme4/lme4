# This code is used to create the testdata called
# test-covariance_structures_data.rds

library(glmmTMB)
library(nlme)
library(lme4)

fm1.glmmTMB <- glmmTMB(Reaction ~ Days + (Days | Subject), sleepstudy, 
                       REML = FALSE)
fm2.glmmTMB <- glmmTMB(Reaction ~ Days + us(Days | Subject), sleepstudy, 
                       REML = FALSE)
fm3.glmmTMB <- glmmTMB(Reaction ~ Days + cs(Days | Subject), sleepstudy, 
                       REML = FALSE)
fm4.glmmTMB <- glmmTMB(Reaction ~ Days + diag(Days | Subject), sleepstudy, 
                       REML = FALSE)

fm1.nlme <- lme(fixed = Reaction ~ Days, random = ~ Days | Subject,
  data = sleepstudy, method = "ML")

fm1.nlme.REML <- lme(fixed = Reaction ~ Days, random = ~ Days | Subject,
                     data = sleepstudy)

fm2.nlme <- lme(fixed = Reaction ~ Days, random = ~ Days | Subject,
                correlation = corCompSymm(form = ~ Days | Subject),
                data = sleepstudy, method = "ML")

fm2.nlme.REML <- lme(fixed = Reaction ~ Days, random = ~ Days | Subject,
                correlation = corCompSymm(form = ~ Days | Subject),
                data = sleepstudy)

# TODO: CHANGE NAME CONVENTION (inconsistent for glmmTMB and nlme...)

# fm1.glmmTMB.diag._vcov -> fm1.glmmTMB.cs_vcov

glmmTMB_mod <- list(
  ## sigma
  fm1.glmmTMB_sigma = sigma(fm1.glmmTMB),
  fm1.glmmTMB.cs_sigma = sigma(fm3.glmmTMB),
  fm1.glmmTMB.diag_sigma = sigma(fm4.glmmTMB),
  fm1.nlme_sigma = sigma(fm1.nlme),
  fm1.nlme.cs_sigma = sigma(fm2.nlme),
  fm1.nlme.REML_sigma = sigma(fm1.nlme.REML),
  fm1.nlme.cs.REML_sigma = sigma(fm2.nlme.REML),
  ## vcov
  fm1.glmmTMB_vcov = vcov(fm1.glmmTMB)$cond,
  fm1.glmmTMB.cs_vcov = vcov(fm3.glmmTMB)$cond,
  fm1.glmmTMB.diag_vcov = vcov(fm4.glmmTMB)$cond,
  fm1.nlme_vcov = vcov(fm1.nlme),
  fm1.nlme.cs_vcov = vcov(fm2.nlme),
  fm1.nlme.REML_vcov = vcov(fm1.nlme.REML),
  fm1.nlme.cs.REML_vcov = vcov(fm2.nlme.REML),
  ## VarCorr
  fm1.glmmTMB_var = VarCorr(fm1.glmmTMB)$cond[[1]],
  fm1.glmmTMB.us_var = VarCorr(fm2.glmmTMB)$cond[[1]],
  fm1.glmmTMB.cs_var = VarCorr(fm3.glmmTMB)$cond[[1]],
  fm1.glmmTMB.diag_var = VarCorr(fm4.glmmTMB)$cond[[1]],
  fm1.nlme_var = getVarCov(fm1.nlme),
  fm1.nlme.REML_var = getVarCov(fm1.nlme.REML),
  fm1.nlme.cs_var = getVarCov(fm2.nlme),
  fm1.nlme.cs.REML_var = getVarCov(fm2.nlme.REML),
  ## ranef
  fm1.glmmTMB_ranef = ranef(fm1.glmmTMB)$cond$Subject,
  fm1.glmmTMB.cs_ranef = ranef(fm3.glmmTMB)$cond$Subject,
  fm1.glmmTMB.diag_ranef = ranef(fm4.glmmTMB)$cond$Subject,
  fm1.nlme_ranef = ranef(fm1.nlme),
  fm1.nlme.REML_ranef = ranef(fm1.nlme.REML),
  fm1.nlme.cs_ranef = ranef(fm2.nlme),
  fm1.nlme.cs.REML_ranef = ranef(fm2.nlme.REML),
  ## predict
  fm1.glmmTMB_predict = predict(fm1.glmmTMB),
  fm1.glmmTMB.cs_predict = predict(fm3.glmmTMB),
  fm1.glmmTMB.diag_predict = predict(fm4.glmmTMB),
  fm1.nlme_predict = predict(fm1.nlme),
  fm1.nlme.REML_predict = predict(fm1.nlme.REML),
  fm1.nlme.cs_predict = predict(fm2.nlme),
  fm1.nlme.cs.REML_predict = predict(fm2.nlme.REML)
)

saveRDS(glmmTMB_mod, "inst/testdata/test-covariance_structures_data.rds")

# What is going on here???

set.seed(1)

all.equal(predict(fm1.glmmTMB), predict(fm1), check.attributes = FALSE)
