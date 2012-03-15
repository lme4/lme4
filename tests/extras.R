## library(lme4)
## This example takes too long
## if (isTRUE(try(data(star, package = 'mlmRev')) == 'star')) {
##     fm1 <- lme4:::carryOver(math ~ gr+sx*eth+cltype+(yrs|id)+(1|tch)+(yrs|sch),
##                             star, yrs ~ tch/id,
##                             control = list(msV = 1, nit = 0, grad = 0))
##     print(fm1, corr = FALSE)
## }
