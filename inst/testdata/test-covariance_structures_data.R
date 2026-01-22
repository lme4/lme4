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

startvec <- c(Asym = 200, xmid = 725, scal = 350)
nm.nlme <- nlme(
  circumference ~ SSlogis(age, Asym, xmid, scal),
  data = Orange,
  fixed = Asym + xmid + scal ~ 1,   
  random = Asym ~ 1 | Tree, 
  start = startvec)

nm.nlme.cs <- nlme(
  circumference ~ SSlogis(age, Asym, xmid, scal),
  data = Orange,
  fixed = Asym + xmid + scal ~ 1,   
  random = Asym ~ 1 | Tree, 
  correlation = corCompSymm(form = Asym ~ 1 | Tree),
  start = startvec)

# tests regarding stdev and corr of glmms

simfun_gamma <- function(ngrp = 50, nrep = 50, shape_gam = 2, intercept = 1, 
                         theta_val = 2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  dd <- expand.grid(group = 1:ngrp, rep = 1:nrep)
  dd$y <- simulate(~ 1 + (1 | group), newdata = dd, 
                   family = Gamma(link = "log"), 
                   newparams = list(
                     theta = theta_val, beta = 1, 
                     sigma = 1/sqrt(shape_gam)))[[1]]
  dd
}

simfun_pois <- function(ngrp = 50, nrep = 50, intercept = 1,
                        theta_val = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  dd <- expand.grid(group = 1:ngrp, rep = 1:nrep)
  
  dd$y <- simulate(~ 1 + (1 | group), newdata = dd,
                   family = poisson(link = "log"), 
                   newparams = list(theta = theta_val, beta = intercept))[[1]]
  dd
}

dd2 <- simfun_gamma(seed = 101)

dd3 <- simfun_pois(seed = 101)

TMB1 <- glmmTMB(y ~ 1 + (1|group), family = Gamma(link = "log"), data = dd2)
TMB_v1 <- VarCorr(TMB1)$cond$group

TMB2 <- glmmTMB(y ~ 1 + (1|group), family = poisson(link = "log"), data = dd3)
TMB_v2 <- VarCorr(TMB2)$cond$group

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
  fm1.glmmTMB_logLik = logLik(fm1.glmmTMB),
  fm1.glmmTMB.cs_logLik = logLik(fm3.glmmTMB),
  fm1.glmmTMB.diag_logLik = logLik(fm4.glmmTMB),
  fm1.glmmTMB.ar1_logLik = logLik(fm5.glmmTMB),
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
  gm.glmmTMB_logLik = logLik(gm.glmmTMB.us),
  gm.glmmTMB.cs_logLik = logLik(gm.glmmTMB.cs),
  gm.glmmTMB.diag_logLik = logLik(gm.glmmTMB.diag),
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
  gm.glmmTMB.diag_predict = predict(gm.glmmTMB.diag),
  ## Third: tests for nlmer
  ## sigma
  nm.nlme_sigma = sigma(nm.nlme),
  nm.nlme.cs_sigma = sigma(nm.nlme.cs),
  ## log likelihoods
  nm.nlme_logLik = logLik(nm.nlme),
  nm.nlme.cs_logLik = logLik(nm.nlme.cs),
  ## vcov
  nm.nlme_vcov = vcov(nm.nlme),
  nm.nlme.cs_vcov = vcov(nm.nlme.cs),
  ## ranef
  nm.nlme_ranef = ranef(nm.nlme),
  nm.nlme.cs_ranef = ranef(nm.nlme.cs),
  ## predict
  nm.nlme_predict = predict(nm.nlme),
  nm.nlme.cs_predict = predict(nm.nlme.cs),
  ## glmm stdevs
  TMB_v1 = TMB_v1,
  TMB_v2 = TMB_v2
)

saveRDS(cov_lmer_test, "inst/testdata/test-covariance_structures_data.rds")
saveRDS(mlmRev::Contraception, "inst/testdata/Contraception.rds")

