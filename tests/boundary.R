## In both of these cases boundary fit (i.e. estimate of zero RE
## variance) is *incorrect*. (Nelder_Mead, restart_edge=FALSE) is the
## only case where we get stuck; either optimizer=bobyqa or
## restart_edge=TRUE (default) works

library(lme4)

## Stephane Laurent:
dat <- read.csv(system.file("testdata","dat20101314.csv",package="lme4"))

fit <- lmer(y ~ (1|Operator)+(1|Part)+(1|Part:Operator), data=dat)
fit_b <- lmer(y ~ (1|Operator)+(1|Part)+(1|Part:Operator), data=dat,
              control=lmerControl(optimizer="bobyqa",restart_edge=FALSE))
fit_c <- lmer(y ~ (1|Operator)+(1|Part)+(1|Part:Operator), data=dat,
              control=lmerControl(restart_edge=FALSE))
## tol=1e-5 seems OK in interactive use but not in R CMD check ... ??
stopifnot(all.equal(getME(fit,"theta"),getME(fit_b,"theta"),tol=2e-5))
stopifnot(all(getME(fit,"theta")>0))

## Manuel Koller

source(system.file("testdata","koller-data.R",package="lme4"))
ldata <- getData(13)
## old (backward compatible/buggy)
fm4 <- lmer(y ~ (1|Var2), ldata, control=lmerControl(use.last.params=TRUE))
fm4b <- lmer(y ~ (1|Var2), ldata, control=lmerControl(restart_edge=FALSE,
                                  use.last.params=TRUE))
stopifnot(getME(fm4b,"theta")==0)
fm4c <- lmer(y ~ (1|Var2), ldata, control=lmerControl(optimizer="bobyqa",
                                  use.last.params=TRUE))
stopifnot(all.equal(getME(fm4,"theta"),getME(fm4c,"theta"),tol=1e-4))
stopifnot(all(getME(fm4,"theta")>0))

## new: doesn't get stuck at edge any more,  but gets stuck somewhere else ...
fm5 <- lmer(y ~ (1|Var2), ldata)
fm5b <- lmer(y ~ (1|Var2), ldata, control=lmerControl(restart_edge=FALSE))
fm5c <- lmer(y ~ (1|Var2), ldata, control=lmerControl(optimizer="bobyqa"))
stopifnot(all(getME(fm4,"theta")>0))

if (FALSE) {
    ## additional stuff for diagnosing Nelder-Mead problems.
    ## library(nloptr) call commented out to avoid R CMD check problems/needing to
    ##  Suggest: nloptr
    
    library(optimx)
    fm5d <- update(fm5,control=lmerControl(optimizer="optimx",
                       optCtrl=list(method="L-BFGS-B")))

    ## library(nloptr)
    defaultControl <- list(algorithm="NLOPT_LN_BOBYQA",xtol_rel=1e-6,maxeval=1e5)
    nloptwrap2 <- function(fn,par,lower,upper,control=list(),...) {
        for (n in names(defaultControl)) 
            if (is.null(control[[n]])) control[[n]] <- defaultControl[[n]]
        res <- nloptr(x0=par,eval_f=fn,lb=lower,ub=upper,opts=control,...)
        with(res,list(par=solution,
                      fval=objective,
                      conv=if (status>0) 0 else status,
                      message=message))
    }
    fm5e <- update(fm5,control=lmerControl(optimizer="nloptwrap2"))

    mList <- setNames(list(fm4,fm4b,fm4c,fm5,fm5b,fm5c,fm5d,fm5e),
                      c("NM/uselast","NM/uselast/norestart","bobyqa/uselast",
                        "NM","NM/norestart","bobyqa","LBFGSB","nloptr/bobyqa"))
    pp <- profile(fm5c,which=1)
    dd <- as.data.frame(pp)
    par(las=1,bty="l")
    v <- sapply(mList,
                function(x) sqrt(VarCorr(x)[[1]]))
    (v2 <- sapply(mList,getME,"theta"))

    plot(.zeta^2~.sig01,data=dd,type="b")
    abline(v=v)

}
