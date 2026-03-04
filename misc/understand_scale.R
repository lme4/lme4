library(lme4)
simfun_gam <- function(ngrp = 50, nrep = 50, shape_gam = 2, intercept = 1, 
                       seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  dd <- expand.grid(group = 1:ngrp, rep = 1:nrep)
  dd$y <- simulate(~ 1 + (1 | group), newdata = dd, 
                     family = Gamma(link = "log"), 
                     newparams = list(
                      theta = 1, beta = 1, sigma = 1/sqrt(shape_gam)))[[1]]
  return(dd) 
}

dd1 <- simfun_gam(seed = 101)

m1 <- glmer(y ~ 1 + (1|group), family = Gamma(link = "log"), data = dd1)

sigma(m1)

VarCorr(m1)

sqrt(deviance(m1) / df.residual(m1))

deviance(m1)
identical(c(-2*logLik(m1)), m1@devcomp$cmp[["dev"]])
identical(lme4:::devCrit(m1), m1@devcomp$cmp[["dev"]])
deviance(m1)
-2*logLik(m1)
## what's the difference between 'sigmaML' and the deviance criterion/obj function?
n <- nrow(dd1)
with(as.list(m1@devcomp$cmp), c(dev, pwrss, sqrt(pwrss/n), sigmaML))

## what is pwrss and why is it different from opt$fval?
## sigmaML is pwrss/n
m1@devcomp$cmp


m1@resp$aic() / df.residual(m1)
-2*logLik(m1) / df.residual(m1)
