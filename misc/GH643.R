library(lme4)
library(glmmTMB)
library(GLMMadaptive)

## what is going on here?  Try with simulated data?
data(epil2)
epil2_red <- subset(epil2, y>0)
form <- y ~ trt + (1|subject)
fam <- Gamma(link = "log")
mod_glmer <- glmer(form, data=epil2_red, family=fam)
mod_glmmTMB <- glmmTMB(form, data=epil2_red, family=fam)
mod_GA <- mixed_model(y ~ trt,
                      random = ~ 1|subject, data=epil2_red,
                      family=fam, control = list(nAGQ=1))

AIC(mod_glmer, mod_glmmTMB, mod_GA)   # very different AIC
## Warning message (ignorable):
## In AIC.default(mod_glmer, mod_glmmTMB, mod_GA) :
##   models are not all fitted to the same number of observations
## why?? nobs() is the same for all three ...
## mixed_model lists nobs as 58 ... == length(unique(epil2$subject) - 1)

## VarCorr ???
## Gamma should report VarCorr as the same as theta[[1]] (unscaled)?
##  double-check 'fix'
mod_list <- list(glmer=mod_glmer, glmmTMB = mod_glmmTMB, GA=mod_GA)
sapply(mod_list, function(x) unlist(fixef(x)[1:2]))

## very different -- if we believe these then
## glmer >> glmmTMB >> GA
sapply(mod_list, function(x) -1*c(logLik(x)))

rssq <- function(x) sum(residuals(x, type = "pearson")^2)
sapply(mod_list,  rssq)
       sapply(mod_list, function(x) unlist(VarCorr(x)))

fixef(mod_glmer)
fixef(mod_glmmTMB)
logLik(mod_glmer)
logLik(mod_glmmTMB)
logLik(mod_glmer)
logLik(mod_GA)

predict_fun <- function(x) {
  if (inherits(x, "MixMod")) {
    return(GLMMadaptive:::predict.MixMod(x, type_pred = "link",
                                  type= "subject_specific"))
  }
  predict(x)
}

## I think we want type = "subject_specific" here?
## predict(mod_GA, type_pred = "link", type = "subject_specific")
## Error in Z * EBs$post_modes[id, seq_len(ncz), drop = FALSE] : 
##   non-conformable arrays

predict_fun(mod_GA)
predict(mod_GA)
pp <- sapply(mod_list, predict_fun)
pairs(pp, gap = 0)

## test likelihood functions directly (with appropriately transformed
## parameters, e.g. theta[glmmTMB] = log(theta[glmer]) ?
lmer_pars <- getME(mod_glmer, c("theta", "beta"))
glmmTMB_pars <- with(mod_glmmTMB$obj$env,
                     last.par.best[-random])

ffun <- function(x) factor(rep(NA, length(x)))
mod_glmmTMB2 <- update(mod_glmmTMB,
       start = list(beta = lmer_pars$beta,
                    theta = log(lmer_pars$theta),
                    betadisp = 0),
       map = list(beta = ffun(lmer_pars$beta),
                  theta = ffun(lmer_pars$theta)))

logLik(mod_glmmTMB2)
logLik(mod_glmmTMB)
logLik(mod_glmer)
sigma(mod_glmmTMB)
sigma(mod_GA)
sigma(mod_glmer)


library(glmmTMB)
library(lme4)

set.seed(123)
N <- 1000
data <- data.frame(y=rgamma(N,3,30),
                   grp=sample(1:20,N,replace=T))
                   
glmmTMB_mod <-glmmTMB(y~1+
                        (1|grp),
                      family=Gamma(link="log"),
                      data=data)

glmer_mod <- glmer(y~1+
                     (1|grp),
                   family=Gamma(link="log"),
                   data=data)

glmer_mod2 <- update(glmer_mod,
                     start = list(beta = fixef(glmmTMB_mod)$cond,
                                  theta = exp(getME(glmmTMB_mod, "theta"))))
glm_mod <- glm(y~1,
               family=Gamma(link="log"),
               data=data)

fixef(glmer_mod)
fixef(glmmTMB_mod)
VarCorr(glmer_mod)
VarCorr(glmmTMB_mod)
## VarCorr are somewhat different (glmmTMB finds a singular fit,
## glmer_mod doesn't

logLik(glmer_mod)
logLik(glmmTMB_mod)

pp <- profile(glmer_mod, which = "theta_") ## can't profile, ugh.
## do it the hard way?

par(mfrow=c(2,2),mar=c(3,3,2,2))
hist(as.matrix(simulate(glmmTMB_mod, nsim = 1))[,1],
     main="glmmTMB",xlab="",breaks=seq(0,1.5,0.05))
hist(as.matrix(simulate(glmer_mod, nsim = 1))[,1],
     main="glmer",xlab="",breaks=seq(0,1.5,0.05))
hist(as.matrix(simulate(glm_mod, nsim = 1))[,1],
     main="glm",xlab="",breaks=seq(0,1.5,0.05))
hist(data$y,main="original distribution",xlab="",breaks=seq(0,1.5,0.05))
