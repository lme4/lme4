## how do we want mkVarCorr to behave for structured matrices, including
## those with big cov matrices?

library(lme4)
library(glmmTMB)

library(lme4)
set.seed(20260217)
dd <- expand.grid(f = factor(1:100), g = factor(1:30))
dd$y <- glmmTMB::simulate_new(~ 1 + ar1(0 + f | g),
                      newdata = dd,
                      newparams = list(beta = 0, theta = c(0,0)),
                      family = gaussian)[[1]]

system.time(
  m <- lmer(y ~ 1 + ar1(0 + f | g), data = dd)
)

##   user  system elapsed 
## 78.425 152.311  15.291 

## note covmat is still full size ...
object.size(VarCorr(m))
str(VarCorr(m))

m <- lmer(y ~ 1 + ar1(0 + f | g), data = dd)
VarCorr(m)
VarCorr(m, full = TRUE)
VarCorr(m, full_cor = FALSE)
undebug(lme4:::VarCorr.merMod)
debug(mkVarCorr)
