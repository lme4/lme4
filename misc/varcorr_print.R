## how do we want mkVarCorr to behave for structured matrices, including
## those with big cov matrices?

library(lme4)
library(glmmTMB)

library(lme4)
set.seed(20260217)
dd <- expand.grid(f = factor(1:30), g = factor(1:30))
dd$y <- glmmTMB::simulate_new(~ 1 + ar1(0 + f | g),
                      newdata = dd,
                      newparams = list(beta = 0, theta = c(0,0)),
                      family = gaussian)[[1]]
        

