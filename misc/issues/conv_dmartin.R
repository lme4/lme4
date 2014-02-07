## From http://stackoverflow.com/questions/21344555/convergence-error-for-development-version-of-lme4
## **tl;dr** this looks like a false positive -- I don't see any
## particularly important differences among the fits with a variety of
## different optimizers, although it looks as though the outliers are
## the built-in Nelder-Mead optimizer and nlminb; built-in bobyqa, and
## bobyqa and Nelder-Mead from the nloptr package, give extremely
## close answers, and no warnings.
## 

## I put the `dput` output in a separate file:

source("convdat.R")

## Run the whole gamut of possible optimizers: built-in N-M and
## bobyqa; nlminb and L-BFGS-B from base R, via the `optimx` package;
## and the `nloptr` versions of N-M and bobyqa.

library(lme4)
g0.bobyqa <- glmer(resp ~ months.c * similarity * percSem +
                 (similarity | subj),
      family = binomial, data = myData,
                   control=glmerControl(optimizer="bobyqa"))
g0.NM <- update(g0.bobyqa,control=glmerControl(optimizer="Nelder_Mead"))
library(optimx)
g0.nlminb <- update(g0.bobyqa,control=glmerControl(optimizer="optimx",
                              optCtrl=list(method="nlminb")))
g0.LBFGSB <- update(g0.bobyqa,control=glmerControl(optimizer="optimx",
                              optCtrl=list(method="L-BFGS-B")))

library(nloptr)
## from https://github.com/lme4/lme4/issues/98:
defaultControl <- list(algorithm="NLOPT_LN_BOBYQA",xtol_rel=1e-6,maxeval=1e5)
nloptwrap2 <- function(fn,par,lower,upper,control=list(),...) {
    for (n in names(defaultControl)) 
      if (is.null(control[[n]])) control[[n]] <- defaultControl[[n]]
    res <- nloptr(x0=par,eval_f=fn,lb=lower,ub=upper,opts=control,...)
    with(res,list(par=solution,
                  fval=objective,
                  feval=iterations,
                  conv=if (status>0) 0 else status,
                  message=message))
}
g0.bobyqa2 <- update(g0.bobyqa,control=glmerControl(optimizer=nloptwrap2))
g0.NM2 <- update(g0.bobyqa,control=glmerControl(optimizer=nloptwrap2,
                           optCtrl=list(algorithm="NLOPT_LN_NELDERMEAD")))

## Summarize results.  We get warnings from `nlminb`, `L-BFGS-B`, and
## Nelder-Mead (but the size of the max abs gradient is largest from
## Nelder-Mead)

getpar <- function(x) c(getME(x,c("theta")),fixef(x))
modList <- list(bobyqa=g0.bobyqa,NM=g0.NM,nlminb=g0.nlminb,
                bobyqa2=g0.bobyqa2,NM2=g0.NM2,LBFGSB=g0.LBFGSB)
ctab <- sapply(modList,getpar)
library(reshape2)
mtab <- melt(ctab)
library(ggplot2)
theme_set(theme_bw())
ggplot(mtab,aes(x=Var2,y=value,colour=Var2))+
    geom_point()+facet_wrap(~Var1,scale="free")

## Just the 'good' fits:

ggplot(subset(mtab,Var2 %in% c("NM2","bobyqa","bobyqa2")),
       aes(x=Var2,y=value,colour=Var2))+
    geom_point()+facet_wrap(~Var1,scale="free")

## Coefficient of variation of estimates among optimizers:

summary(cvvec <- apply(ctab,1,function(x) sd(x)/mean(x)))

## The highest CV is for `months.c`, which is still only about 4% ...

## The log-likelihoods don't differ very much: NM2 gives the max
## log-likelihood, and all the 'good' ones are very close (even the
## 'bad' ones are at most 1% different)

likList <- sapply(modList,logLik)
round(log10(max(likList)-likList),1)
##  bobyqa      NM  nlminb bobyqa2     NM2  LBFGSB 
##    -8.5    -2.9    -2.0   -11.4    -Inf    -5.0 

