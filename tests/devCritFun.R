library(lme4)

## ----------------------------------------------------------------------
## test that deviance(REMLfit, REML = FALSE) gives the same answer as
## the ML objective function at the REML fit
## ----------------------------------------------------------------------
set.seed(1)
w <- runif(nrow(sleepstudy))
fm <- lmer(Reaction ~ Days + (Days | Subject),
           sleepstudy, weights = w)
dfun <- update(fm, devFunOnly = TRUE, REML = FALSE)
stopifnot(all.equal(deviance(fm, REML = FALSE),
                    dfun(getME(fm, "theta"))))

## ----------------------------------------------------------------------
## TODO: test the opposite case that deviance(MLfit, REML = TRUE)
## gives the same answer as the REML objective function at the ML fit
## ----------------------------------------------------------------------
