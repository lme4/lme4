library(lme4)
dd  <- data.frame(g = factor(rep(1:3, each = 3)))
dd$y <- simulate(~ 1 + (1|g),
              family = poisson,
              seed = 101,
              newdata = dd,
              newparams = list(beta = 0, theta = 2))[[1]]

## constructive::construct(dd)
## equivalent, maybe helpful in porting to Julia etc.
dd <- data.frame(
  g = factor(rep(1:3, each = 3L)),
  y = as.integer(c(0, 0, 1, 3, 5, 4, 0, 1, 0))
)

g1 <- glmer(y ~ 1 + (1|g), data = dd, family = poisson, control = glmerControl(nAGQ0initStep = FALSE))
g2 <- update(g1, control = glmerControl(compDev = FALSE, nAGQ0initStep = FALSE))

stopifnot(identical(VarCorr(g1), VarCorr(g2)))

dfun1 <- getME(g1, "devfun")
dfun2 <- getME(g2, "devfun")
stopifnot(!identical(environment(dfun1), environment(dfun2)))

pars <- c(1.16204380240176, -0.183842216759865)
stopifnot(identical(dfun1(pars), dfun2(pars)))

stopifnot(!get("compDev", environment(dfun2)))

## compare with lme4pureR (not identical)
remotes::install_github("lme4/lme4pureR")
library(lme4pureR)

ll <- plsform(formula(g2), data = dd, family = poisson)
devf <- do.call(pirls, c(ll, list(family=poisson)))
devf(pars)
all.equal(devf(pars), dfun1(pars), tolerance = 0)
## differs by 6e-6 (could compare in more detail)

## debug(dfun2)
## dfun2(pars)
## ... debug(pwrssUpdate) ## in environment(dfun2)
debug(environment(dfun2)$pwrssUpdate)
debug(lme4:::RglmerWrkIter)
dfun2(pars)
