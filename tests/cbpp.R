library(lme4)
library(nloptwrap)
library(lme4pureR)
form <- cbind(incidence, size - incidence) ~ period + (1 | herd)
glmod <- glFormula(form, cbpp, binomial)
gm1 <- glm(nobars(form), binomial, cbpp)
devf <- pirls(glmod, gm1$y, binomial, weights=gm1$prior.weights,
              eta=gm1$linear.predictors, nAGQ=1L)
cc <- coef(gm1)
nlminb(c(1,cc), devf, lower=c(0,rep.int(-Inf,length(cc))))[c("par","objective")]
(gm2 <- glmer(form, cbpp, binomial, verbose=2L,start=list(theta=1,fixef=cc)))
dd <- glmer(form, cbpp, binomial, nAGQ=1L, devFunOnly=TRUE)
