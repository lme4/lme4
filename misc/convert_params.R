library(lme4)
library(glmmTMB)

fit_cs_glmm <- glmmTMB(Reaction ~ Days + cs(Days | Subject), sleepstudy, REML = FALSE)
fit_cs <- lmer(Reaction ~ Days + cs(Days | Subject), sleepstudy, REML = FALSE)


getME(fit_cs, "theta")
getME(fit_cs_glmm, "theta")

debug(lme4:::VarCorr.merMod)
VarCorr(fit_cs)
        
lme4:::cs_theta_to_rho
lme4:::cs_rho_to_theta

getME(fit_cs, "theta")
glmmTMB_theta <- getME(fit_cs_glmm, "theta")
m <- matrix(1, 2, 2)
m[1, 2] <- m[2,1] <- lme4:::cs_theta_to_rho(glmmTMB_theta[3], 2)
sdvec <- exp(glmmTMB_theta[1:2])
sqrt(diag(outer(sdvec,  sdvec) * m))
## [1] 0.9 2618078 0.22337456 0.07557558


glmmTMB_to_lme4_cs <- function(fit) {
  theta <- getME(fit, "theta")
  s <- sigma(fit)
  c(exp(theta[1:2])/s, theta[3])
}
glmmTMB_to_lme4_cs(fit_cs_glmm)
