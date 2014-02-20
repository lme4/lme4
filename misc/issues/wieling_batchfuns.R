nloptwrap2 <- function(fn,par,lower,upper,control=list(),...) {
    ## use DMB's suggested setting of xtol_rel=1e-6
    if (is.null(control$xtol_rel)) control$xtol_rel <- 1e-6
    if (is.null(control$maxeval)) control$maxeval <- 1e5
    res <- nloptr(x0=par,eval_f=fn,lb=lower,ub=upper,opts=control,...)
    with(res,list(par=solution,
                  fval=objective,
                  feval=iterations,
                  conv=if (status>0) 0 else status,
                  message=message))
}

fullForm <- PronDistStdDutch.c ~ Geo + PopSize.log_residGeo.z +
    PopAvgAge_residPopAvgIncome.log_Geo.z +
    PopAvgIncome.log_residGeo.z + WordFreq.log.z +
    WordIsNounOrAdverb + WordVCratio.log.z +
    (1 | Word) +
    (0 + PopSize.log_residGeo.z + PopAvgAge_residPopAvgIncome.log_Geo.z + PopAvgIncome.log_residGeo.z | Word) +
    (1 + WordFreq.log.z + WordIsNounOrAdverb | Location) +
    (1 | Transcriber)

fitLme4.0 <- function(data=subdat) {
    tval <- system.time(fit <- lme4.0::lmer(fullForm, data=data))
    list(time=tval, fit=fit)
}

fitLme4 <- function(optCtrl=list(), optimizer, data=subdat) {
    tval <- system.time(fit <- lme4::lmer(fullForm,data=data,
          control=lmerControl(optimizer=optimizer, optCtrl=optCtrl)))
    list(time=tval, fit=fit)
}

## use this rather than the nloptwrap2 control above so that
## the control settings are recorded in the fitted object ...
nlopt <- function(algorithm,xtol_rel=1e-6,maxeval=1e5) {
    list(algorithm=algorithm,xtol_rel=xtol_rel,maxeval=maxeval)
}
argList <-
    list(
        nm1=list(optimizer="Nelder_Mead"),
        nm2=list(optimizer="Nelder_Mead",optCtrl=list(maxfun=1e5)),
        bobyqa1=list(optimizer="bobyqa"),
        nlminb=list(optimizer="optimx",optCtrl=list(method="nlminb")),
        lbfgsb=list(optimizer="optimx",optCtrl=list(method="L-BFGS-B")),
        ## nloptr derivative-free choices
        nloptbobyqa=list(optimizer="nloptwrap2",
            optCtrl=nlopt("NLOPT_LN_BOBYQA")),
        nloptcobyla=list(optimizer="nloptwrap2",
            optCtrl=nlopt("NLOPT_LN_COBYLA")),
        nloptNM=list(optimizer="nloptwrap2",
            optCtrl=nlopt("NLOPT_LN_NELDERMEAD")),
        nloptsubplex=list(optimizer="nloptwrap2",
            optCtrl=nlopt("NLOPT_LN_SBPLX"))
        )
