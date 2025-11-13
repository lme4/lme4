# This code is used to create the testdata called
# test-covariance_structures_data.rds

library(glmmTMB)
library(nlme)
library(lme4)

## First set of tests: lmer
fm1.glmmTMB <- glmmTMB(Reaction ~ Days + (Days | Subject), sleepstudy, 
                       REML = FALSE)
fm2.glmmTMB <- glmmTMB(Reaction ~ Days + us(Days | Subject), sleepstudy, 
                       REML = FALSE)
fm3.glmmTMB <- glmmTMB(Reaction ~ Days + cs(Days | Subject), sleepstudy, 
                       REML = FALSE)
fm4.glmmTMB <- glmmTMB(Reaction ~ Days + diag(Days | Subject), sleepstudy, 
                       REML = FALSE)
sleepstudy$Daysf <- factor(sleepstudy$Days, ordered = TRUE)
fm5.glmmTMB <- glmmTMB(Reaction ~ Daysf + ar1(0 + Daysf | Subject), 
                       sleepstudy, REML = FALSE)

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

## Second set of tests: glmer
data(Contraception, package = "mlmRev")

gm.glmmTMB.us <- glmmTMB(use ~ age + urban + us(1 + urban | district),
                         data = Contraception,
                         family = binomial)
gm.glmmTMB.cs <- glmmTMB(use ~ age + urban + cs(1 + urban | district),
                         data = Contraception,
                         family = binomial)
gm.glmmTMB.diag <- glmmTMB(use ~ age + urban + diag(1 + urban | district),
                           data = Contraception,
                           family = binomial)

# fm1.glmmTMB.diag._vcov -> fm1.glmmTMB.cs_vcov

cov_lmer_test <- list(
  ## First: set of tests for lmer
  ## sigma
  fm1.glmmTMB_sigma = sigma(fm1.glmmTMB),
  fm1.glmmTMB.cs_sigma = sigma(fm3.glmmTMB),
  fm1.glmmTMB.diag_sigma = sigma(fm4.glmmTMB),
  fm1.glmmTMB.ar1_sigma = sigma(fm5.glmmTMB),
  fm1.nlme_sigma = sigma(fm1.nlme),
  fm1.nlme.cs_sigma = sigma(fm2.nlme),
  fm1.nlme.REML_sigma = sigma(fm1.nlme.REML),
  fm1.nlme.cs.REML_sigma = sigma(fm2.nlme.REML),
  ## log likelihoods
  fm1.glmmTMB_logLik = -fm1.glmmTMB$fit$objective,
  fm1.glmmTMB.cs_logLik = -fm3.glmmTMB$fit$objective,
  fm1.glmmTMB.diag_logLik = -fm4.glmmTMB$fit$objective,
  fm1.glmmTMB.ar1_logLik = -fm5.glmmTMB$fit$objective,
  fm1.nlme_logLik = logLik(fm1.nlme),
  fm1.nlme.cs_logLik = logLik(fm2.nlme),
  fm1.nlme.REML_logLik = logLik(fm1.nlme.REML),
  fm1.nlme.cs.REML_logLik = logLik(fm2.nlme.REML),
  ## vcov
  fm1.glmmTMB_vcov = vcov(fm1.glmmTMB)$cond,
  fm1.glmmTMB.cs_vcov = vcov(fm3.glmmTMB)$cond,
  fm1.glmmTMB.diag_vcov = vcov(fm4.glmmTMB)$cond,
  fm1.glmmTMB.ar1_vcov = vcov(fm5.glmmTMB)$cond,
  fm1.nlme_vcov = vcov(fm1.nlme),
  fm1.nlme.cs_vcov = vcov(fm2.nlme),
  fm1.nlme.REML_vcov = vcov(fm1.nlme.REML),
  fm1.nlme.cs.REML_vcov = vcov(fm2.nlme.REML),
  ## VarCorr
  fm1.glmmTMB_var = VarCorr(fm1.glmmTMB)$cond[[1]],
  fm1.glmmTMB.us_var = VarCorr(fm2.glmmTMB)$cond[[1]],
  fm1.glmmTMB.cs_var = VarCorr(fm3.glmmTMB)$cond[[1]],
  fm1.glmmTMB.diag_var = VarCorr(fm4.glmmTMB)$cond[[1]],
  fm1.glmmTMB.ar1_var = VarCorr(fm5.glmmTMB)$cond[[1]],
  fm1.nlme_var = getVarCov(fm1.nlme),
  fm1.nlme.REML_var = getVarCov(fm1.nlme.REML),
  fm1.nlme.cs_var = getVarCov(fm2.nlme),
  fm1.nlme.cs.REML_var = getVarCov(fm2.nlme.REML),
  ## ranef
  fm1.glmmTMB_ranef = ranef(fm1.glmmTMB)$cond$Subject,
  fm1.glmmTMB.cs_ranef = ranef(fm3.glmmTMB)$cond$Subject,
  fm1.glmmTMB.diag_ranef = ranef(fm4.glmmTMB)$cond$Subject,
  fm1.glmmTMB.ar1_ranef = ranef(fm5.glmmTMB)$cond$Subject,
  fm1.nlme_ranef = ranef(fm1.nlme),
  fm1.nlme.REML_ranef = ranef(fm1.nlme.REML),
  fm1.nlme.cs_ranef = ranef(fm2.nlme),
  fm1.nlme.cs.REML_ranef = ranef(fm2.nlme.REML),
  ## predict
  fm1.glmmTMB_predict = predict(fm1.glmmTMB),
  fm1.glmmTMB.cs_predict = predict(fm3.glmmTMB),
  fm1.glmmTMB.diag_predict = predict(fm4.glmmTMB),
  fm1.glmmTMB.ar1_predict = predict(fm5.glmmTMB),
  fm1.nlme_predict = predict(fm1.nlme),
  fm1.nlme.REML_predict = predict(fm1.nlme.REML),
  fm1.nlme.cs_predict = predict(fm2.nlme),
  fm1.nlme.cs.REML_predict = predict(fm2.nlme.REML),
  ## Second: tests for glmer
  ## sigma
  gm.glmmTMB_sigma = sigma(gm.glmmTMB.us),
  gm.glmmTMB.cs_sigma = sigma(gm.glmmTMB.cs),
  gm.glmmTMB.diag_sigma = sigma(gm.glmmTMB.diag),
  ## log likelihoods
  gm.glmmTMB_logLik = -gm.glmmTMB.us$fit$objective,
  gm.glmmTMB.cs_logLik = -gm.glmmTMB.us$fit$objective,
  gm.glmmTMB.diag_logLik = -gm.glmmTMB.diag$fit$objective,
  ## vcov
  gm.glmmTMB_vcov = vcov(gm.glmmTMB.us)$cond,
  gm.glmmTMB.cs_vcov = vcov(gm.glmmTMB.cs)$cond,
  gm.glmmTMB.diag_vcov = vcov(gm.glmmTMB.diag)$cond,
  ## VarCorr
  gm.glmmTMB_var = VarCorr(gm.glmmTMB.us)$cond[[1]],
  gm.glmmTMB.cs_var = VarCorr(gm.glmmTMB.cs)$cond[[1]],
  gm.glmmTMB.diag_var = VarCorr(gm.glmmTMB.diag)$cond[[1]],
  ## ranef
  gm.glmmTMB_ranef = ranef(gm.glmmTMB.us)$cond$Subject,
  gm.glmmTMB.cs_ranef = ranef(gm.glmmTMB.cs)$cond$Subject,
  gm.glmmTMB.diag_ranef = ranef(gm.glmmTMB.diag)$cond$Subject,
  ## predict
  gm.glmmTMB_predict = predict(gm.glmmTMB.us),
  gm.glmmTMB.cs_predict = predict(gm.glmmTMB.cs),
  gm.glmmTMB.diag_predict = predict(gm.glmmTMB.diag)
)

saveRDS(cov_lmer_test, "inst/testdata/test-covariance_structures_data.rds")
saveRDS(mlmRev::Contraception, "inst/testdata/Contraception.rds")

# What is going on here???

set.seed(1)

all.equal(predict(fm1.glmmTMB), predict(fm1), check.attributes = FALSE)
