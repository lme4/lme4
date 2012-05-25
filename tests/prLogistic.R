## data set and formula extracted from ?prLogisticDelta example
##   (Thailand, clustered-data) in prLogistic package
load("prLogistic.RData")
library(lme4)
## FIXME: un-try() after PIRLS failure addressed
try(glmer(rgi~  sex + pped + (1|schoolid), 
                data = dataset, family=binomial))


if (FALSE) {
  library(ggplot2)
  dataset$schoolid <- factor(dataset$schoolid)
  length(levels(dataset$schoolid))
  ttab <- with(dataset,table(schoolid,sex,rgi,pped))
  ttab2 <- with(dataset,table(tapply(rgi,schoolid,
                                     function(x) ifelse(all(x==0),0,
                                                        ifelse(all(x==1),1,0.5)))))
  ggplot(dataset,aes(x=sex,y=rgi,colour=factor(pped)))+stat_sum(alpha=0.5)

  ## library(glmmML)
  ## glmmML_fit <- glmmML(rgi~  sex + pped , cluster=schoolid,
  ##          data = dataset, family=binomial)
  ## glmmML_est <- list(sigma=glmmML_fit$sigma,beta=coef(glmmML_fit))
  ## dput(glmmML_est)
  glmmML_est <- structure(list(sigma = 1.25365353546143,
                               beta = structure(c(-2.19478801858317, 
                               0.548884468743364, -0.623835613907385), .Names = c("(Intercept)", 
                                                                       "sex", "pped"))),
                          .Names = c("sigma", "beta"))

  ## library(lme4.0)
  ## lme4.0_fit <- glmer(rgi~  sex + pped + (1|schoolid), 
  ##                     data = dataset, family=binomial)
  ## lme4.0_est <- list(sigma=unname(sqrt(unlist(VarCorr(lme4.0_fit)))),
  ##                   beta=fixef(lme4.0_fit))
  ##  dput(lme4.0_est)
  ## detach("package:lme4.0")
  lme4.0_est <- structure(list(sigma = 1.25369539060849, beta = structure(c(-2.19474529099587, 
                      0.548900267825802, -0.623934772981894), .Names = c("(Intercept)", 
                          "sex", "pped"))), .Names = c("sigma", "beta"))
  
  devfun <- glmer(rgi~  sex + pped + (1|schoolid), 
                  data = dataset, family=binomial,
                  devFunOnly=TRUE)
  with(glmmML_est,devfun(c(sigma,beta)))
  with(glmmML_est,devfun(c(sigma,beta))) ## FAILS
  devfun <- glmer(rgi~  sex + pped + (1|schoolid), 
                  data = dataset, family=binomial,
                  devFunOnly=TRUE)
  with(lme4.0_est,devfun(c(sigma,beta)))  ## 6326.456
  g3 <- with(glmmML_est,glmer(rgi~  sex + pped + (1|schoolid), 
              data = dataset, family=binomial,
              start=list(theta=sigma,fixef=beta)),
             verbose=10)
}
