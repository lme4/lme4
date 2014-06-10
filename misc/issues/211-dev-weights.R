## @stevencarlislewalker : Good you found the weights  problem!
## ....

## However, the weights problem is still there ... I think release critical
## (Steven accidentally did not use `dfunW`). Here's the updated code. as self-contained R script:
##-----------------------------------------------------------------------

## Now with weights {started from scratch; so self-contained}:
library(lme4)
fm <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
set.seed(1)
w <- runif(nrow(sleepstudy))
fmW <- update(fm, weights = w)
dfunW <- update(fmW, devFunOnly = TRUE, REML = FALSE)
rbind(dev = deviance(fmW, REML = FALSE),
      dfun= dfunW(getME(fmW, "theta")),
      crit= REMLcrit(fmW))
## dev  1612.564
## dfun 1772.279
## crit 1764.115

##-----------------------------------------------------------------------
## where I still wonder why the first two are not the same
## (the third one, `REMLcrit()`,  should differ of course).
