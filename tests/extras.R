library(lme4)
## This example takes  long : only for testLevel >= 3 :
d.ok <- isTRUE(try(data(star, package = 'mlmRev')) == 'star')
if(!interactive() && (lme4:::testLevel() < 3 || !d.ok))
    q("no")

## This worked in an *older* version of lme4.0 
## fm1 <- lme4:::carryOver(math ~ gr+sx*eth+cltype+(yrs|id)+(1|tch)+(yrs|sch),
##                         star, yrs ~ tch/id,
##                         control = list(msV = 1, nit = 0, grad = 0))

system.time(
    fm1 <- lmer(math ~ gr + sx*eth + cltype + schtype + hdeg + clad + exp + trace +
                (yrs | id) + (1 | tch) + (yrs | sch), data = star, verbose = TRUE)
)
##   user  system elapsed 
## 34.991   0.037  35.132 -- lme4.0
## 36.599   0.031  36.745 -- lme4 {bobyqa; 2014-01-09 @ lynne}

sm1 <- summary(fm1)
print(sm1, corr=TRUE, symbolic.cor=TRUE)# now message *and* gives the correlation

