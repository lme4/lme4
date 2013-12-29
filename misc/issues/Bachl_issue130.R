load("Bachl_data.RData")
library(lme4.0)
t0 <- system.time(m0 <- lmer(rtr2 ~ turnsec + (turnsec | kombiid) +
                             (turnsec | turnid) + (turnsec | idnr),
                             verbose = TRUE, d1))
save(list=ls(pattern="^[mt]"),file="Bachl_out.RData")
detach("package:lme4.0")
library("lme4")
t1 <- system.time(m1 <- lmer(rtr2 ~ turnsec + (turnsec | kombiid) +
                             (turnsec | turnid) + (turnsec | idnr),
                             verbose = TRUE, d1))
save(list=ls(pattern="^[mt]"),file="Bachl_out.RData")
t2 <- system.time(m2 <- update(m1,control=lmerControl(optimizer="bobyqa")))
save(list=ls(pattern="^[mt]"),file="Bachl_out.RData")
library(nloptr)
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

t3 <- system.time(m3 <- update(m1,control=lmerControl(optimizer="nloptwrap2")))
save(list=ls(pattern="^[mt]"),file="Bachl_out.RData")
t4 <- system.time(m4 <- update(m1,control=lmerControl(optimizer="nloptwrap2",
                                  optCtrl=list(algorithm="NLOPT_LN_NELDERMEAD"))))
save(list=ls(pattern="^[mt]"),file="Bachl_out.RData")


