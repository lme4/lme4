library(lme4)
library(optimx)
library(lme4.0)
library(nloptr)

## wrapper to make nloptwrap look the way glmer wants it to
nloptwrap2 <- function(fn,par,lower,upper,control,...) {
    res <- nloptr(x0=par,eval_f=fn,lb=lower,ub=upper,opts=control,...)
    with(res,list(par=solution,
                  fval=objective,
                  conv=if (status>0) 0 else status,
                  message=message))
}

set.seed(101)
L <- load("dialectNL.rda")

fullForm <- PronDistStdDutch.c ~ Geo + PopSize.log_residGeo.z +
    PopAvgAge_residPopAvgIncome.log_Geo.z +
    PopAvgIncome.log_residGeo.z + WordFreq.log.z +
    WordIsNounOrAdverb + WordVCratio.log.z +
    (1 | Word) +
    (0 + PopSize.log_residGeo.z +
     PopAvgAge_residPopAvgIncome.log_Geo.z +
     PopAvgIncome.log_residGeo.z | Word) +
    (1 + WordFreq.log.z +
     WordIsNounOrAdverb | Location) +
    (1 | Transcriber)

sampfrac <- 0.01
ntot <- nrow(dialectNL)
subdat <- dialectNL[sample(round(sampfrac*ntot)),]


fitLme4.0 <- function(data=subdat) {
    tval <- system.time(fit <- lme4.0::lmer(fullForm,data=data))
    list(time=tval,fit=fit)
}

fitList0 <- fitLme4.0()

fitLme4 <- function(optCtrl=list(),optimizer,data=subdat) {
    tval <- system.time(fit <- lme4::lmer(fullForm,data=data,
          control=lmerControl(optimizer=optimizer,
                optCtrl=optCtrl)))
    list(time=tval,fit=fit)
}

argList <-
    list(
        nm1=list(optimizer="Nelder_Mead"),
        nm2=list(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e5)),
        bobyqa1=list(optimizer="bobyqa"),
        nlminb=list(optimizer="optimx",optCtrl=list(method="nlminb")),
        lbfgsb=list(optimizer="optimx",optCtrl=list(method="L-BFGS-B")),
        ## nloptr derivative-free choices
        ## (these seem to finish very quickly but I think are mostly way off)
        ## use default xtol_rel=1e-4 throughout to avoid warnings
        nloptbobyqa=list(optimizer="nloptwrap2",
        optCtrl=list(algorithm="NLOPT_LN_BOBYQA",xtol_rel=1e-4)),
        nloptcobyla=list(optimizer="nloptwrap2",
        optCtrl=list(algorithm="NLOPT_LN_COBYLA",xtol_rel=1e-4)),
        nloptNM=list(optimizer="nloptwrap2",
        optCtrl=list(algorithm="NLOPT_LN_NELDERMEAD",xtol_rel=1e-4)),
        nloptsubplex=list(optimizer="nloptwrap2",
        optCtrl=list(algorithm="NLOPT_LN_SBPLX",xtol_rel=1e-4)))

fitList1 <- lapply(argList,
                   function(L) do.call(fitLme4,L))
save(c("fitList0","fitList1"),
     file="wieling_batch1.RData")

