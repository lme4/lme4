L <- load("OddSets1.Rdata")
library(lme4)

## replicate results
selmod_data <- lapply(selmod,"[[","data")
mod_list <- lapply(selmod_data,lmer,
                   formula=dQTcF ~ C + TP + trt -1 +(1 + C|id))
sapply(mod_list,function(x) fixef(x)["C"])
source(system.file("utils", "allFit.R", package="lme4"))
## ugh. Having used lapply(), have to hack the @call slot so update() works ...
## would be easier to try to fix this and/or use a for loop in the first
## place
fit_all_list <- list()
for (i in seq_along(mod_list)) {
    mod_list[[i]]@call[[1]] <- quote(lmer)
    mod_list[[i]]@call$data <- bquote(selmod_data[[.(i)]])
    fit_all_list[[i]] <- allFit(mod_list[[i]])
}
## all estimates are *mostly* similar across optimizers, but not
## entirely ... nloptwrap gives different (smaller) answers
## for 2,4,5,6,9,10
sapply(fit_all_list,
       function(x) summary(x)$fixef[,"C"])
## check conv. warnings: nloptwrap complains a lot but gives correct(er?)
## answers
sapply(fit_all_list,
       function(x) sapply(summary(x)$msgs,length))

## try other tools ... glmmTMB, INLA??

library(glmmTMB)
lapply(selmod_data,glmmTMB,
       formula=dQTcF ~ C + TP + trt -1 +(1 + C|id))
## evaluation problem with lapply()! need to track that down ...
mod_list_2 <- list()
for (i in seq_along(selmod_data)) {
    mod_list_2 <- c(mod_list_2,
                    list(glmmTMB(selmod_data[[i]],
                                 formula=dQTcF ~ C + TP + trt -1 +(1 + C|id))))
}
round(sapply(mod_list_2,function(x) fixef(x)$cond["C"]),3)

library(INLA)
i1 <- inla(formula=dQTcF ~ C + TP + trt -1 +f(id,C,model="iid"),
     data=selmod_data[[1]])
summary(i1)$fixed["C","mean"]

## suspect glmmADMB will be slow & flaky for this ...
if (FALSE) {
    library(glmmADMB)
    lapply(selmod_data,glmmadmb,
           formula=dQTcF ~ C + TP + trt -1 +(1 + C|id),
           family="gaussian")
}

## thoughts/TO DO:
##   compare VarCorr estimates?
##   ditto, BLUPs/conditional modes
##   look at slices -- is this an optimization problem?  (Doesn't
##       seem so, based on allFit results)
##   look at standard errors of estimates too
##  (improve broom to handle all of these cases??)
