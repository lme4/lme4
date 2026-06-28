library(lme4)
dd  <- data.frame(g = factor(rep(1:3, each = 3)))
y <- simulate(~ 1 + (1|g),
              family = poisson,
              seed = 101,
              newdata = dd,
              newparams = list(beta = 0, theta = 2))[[1]]
g1 <- glmer(y ~ 1 + (1|g), data = dd, family = poisson, control = glmerControl(nAGQ0initStep = FALSE))
g2 <- update(g1, control = glmerControl(compDev = FALSE, nAGQ0initStep = FALSE))

stopifnot(identical(VarCorr(g1), VarCorr(g2)))

dfun1 <- getME(g1, "devfun")
dfun2 <- getME(g2, "devfun")
stopifnot(!identical(environment(dfun1), environment(dfun2)))

pars <- c(1.16204380240176, -0.183842216759865)
stopifnot(identical(dfun1(pars), dfun2(pars)))

stopifnot(!get("compDev", environment(dfun2)))
## debug(dfun2)
## dfun2(pars)
## ... debug(pwrssUpdate) ## in environment(dfun2)
debug(environment(dfun2)$pwrssUpdate)
debug(lme4:::RglmerWrkIter)
dfun2(pars)
