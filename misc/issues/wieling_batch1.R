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

( load("dialectNL.rda") ) # "dialectNL
str(dialectNL)
## Find more factors: which are integer or half-integer?
which(isI <- sapply(dialectNL, function(x) !is.factor(x) &&
                    isTRUE(all.equal(2 * x, round(2 * x), tol=1e-12))))
## five of them
## whereas this one, interestingly has quite interesting values:
with(dialectNL, table(SpeakerIsMale))
## SpeakerIsMale
##      0   0.33    0.5   0.75      1 
##  60708    555   3243    513 150082 

dialectNL[isI] <- lapply(dialectNL[isI], as.factor)
str(dialectNL)

fullForm <- PronDistStdDutch.c ~ Geo + PopSize.log_residGeo.z +
    PopAvgAge_residPopAvgIncome.log_Geo.z +
    PopAvgIncome.log_residGeo.z + WordFreq.log.z +
    WordIsNounOrAdverb + WordVCratio.log.z +
    (1 | Word) +
    (0 + PopSize.log_residGeo.z + PopAvgAge_residPopAvgIncome.log_Geo.z + PopAvgIncome.log_residGeo.z | Word) +
    (1 + WordFreq.log.z + WordIsNounOrAdverb | Location) +
    (1 | Transcriber)

sampfrac <- 0.01
ntot <- nrow(dialectNL)
set.seed(101)
## random subset of 1% :
subdat <- dialectNL[sample(round(sampfrac*ntot)),]

fitLme4.0 <- function(data=subdat) {
    tval <- system.time(fit <- lme4.0::lmer(fullForm, data=data))
    list(time=tval, fit=fit)
}

fitList0 <- fitLme4.0()

fitLme4 <- function(optCtrl=list(), optimizer, data=subdat) {
    tval <- system.time(fit <- lme4::lmer(fullForm,data=data,
          control=lmerControl(optimizer=optimizer, optCtrl=optCtrl)))
    list(time=tval, fit=fit)
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

if(FALSE){ ## The only reason to use plyr is the .progress option, 
           ## but that's very inaccurate  
    library(plyr)
    fitList1 <- llply(argList,
                      function(L) do.call(fitLme4,L),.progress="text")
} else
fitList1 <- lapply(argList, do.call, what=fitLme4)
save("fitList0","fitList1",
     file="wieling_batch1.RData")

fitList <- c(lme4.0 = list(fitList0), fitList1)
cat("Timing:\n")
t(sapply(fitList, `[[`, 1)[1:3,])
fits <- lapply(fitList, `[[`, 2)
str(fits, max=1)
ff <- sapply(fits, fixef)
llik <- sapply(fits, logLik)
ff1 <- ff; rownames(ff1) <- abbreviate(rownames(ff))
ff1 <- cbind(logLik=llik, t(ff1))[order(llik, decreasing=TRUE),]
op <- options(width=99, digits=4); ff1 ;  options(op)
##                   logLik      (In)     Geo     PS._     PAA_     PAI.    WF..    WINO    WVC.
## lme4.0            -168.8 -0.009019  1.2892 -0.01601  0.07971 -0.03966 0.01742 0.02257 0.06792
## nlminb.REML       -168.8 -0.008997  1.2555 -0.01391  0.06986 -0.03367 0.01742 0.02257 0.06792
## bobyqa1.REML      -168.8  0.075895 -0.2183  0.27273 -0.51582  0.31628 0.01742 0.02257 0.06792
## lbfgsb.REML       -168.8 -0.009024  1.2893 -0.01600  0.07973 -0.03966 0.01742 0.02256 0.06793
## nm2.REML          -175.6 -0.586885  9.7266 -2.66368  4.58369 -3.74450 0.02264 0.05388 0.06744
## nm1.REML          -175.6 -0.578448  9.7740 -2.62885  4.56213 -3.71126 0.02254 0.05400 0.06745
## nloptsubplex.REML -201.2 -0.006643  1.3084 -0.01446  0.07952 -0.04440 0.02049 0.03881 0.06168
## nloptcobyla.REML  -204.4 -0.002736  1.2456  0.03927  0.01427  0.02277 0.02373 0.05002 0.06419
## nloptNM.REML      -217.8  0.052678  1.2023  0.22124 -0.16660  0.19470 0.02397 0.04735 0.06162
## nloptbobyqa.REML    -Inf  0.013847  0.8769  0.11600 -0.14256  0.14017 0.02163 0.04853 0.06347

## Note that  'bobyqa' (# 3) gets the same good logLik, but 
## quite different coefficients ... but then these seem not to be significant..
