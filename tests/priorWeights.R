rm(list = ls())
set.seed(1)
nGroups <- 100
nObs <- 1000

# explanatory variable with a fixed effect
explVar1 <- rnorm(nObs)
explVar2 <- rnorm(nObs)

# random intercept among levels of a grouping factor
groupFac <- as.factor(rep(1:nGroups,each=nObs/nGroups))
randEff0 <- rep(rnorm(nGroups),each=nObs/nGroups)
randEff1 <- rep(rnorm(nGroups),each=nObs/nGroups)
randEff2 <- rep(rnorm(nGroups),each=nObs/nGroups)

# residuals with heterogeneous variance
residSD <- rpois(nObs,1)+1
#residSD <- rep(1,nObs)
residError <- rnorm(nObs,sd=residSD)

# response variable
respVar <- randEff0 + (1+randEff1)*explVar1 + (1+randEff2)*explVar2 + residError

# rename to fit models on one line
y <- respVar
x <- explVar1
z <- explVar2
g <- groupFac
v <- residSD^2
w <- 1/v

library("nlme")
lmeMods <- list(
    ML1 = lme(y ~ x, random = ~ 1|g, weights = varFixed(~v), method = "ML"),
    REML1 = lme(y ~ x, random = ~ 1|g, weights = varFixed(~v), method = "REML"),
    ML2 = lme(y ~ x, random = ~ x|g, weights = varFixed(~v), method = "ML"),
    REML2 = lme(y ~ x, random = ~ x|g, weights = varFixed(~v), method = "REML"),
    ML1 = lme(y ~ x+z, random = ~ x+z|g, weights = varFixed(~v), method = "ML"),
    REML2 = lme(y ~ x+z, random = ~ x+z|g, weights = varFixed(~v), method = "REML"))

library("lme4")
lmerMods <- list(
    ML1 = lmer(y ~ x + (1|g), weights = w, REML = FALSE),
    REML1 = lmer(y ~ x + (1|g), weights = w, REML = TRUE),
    ML2 = lmer(y ~ x + (x|g), weights = w, REML = FALSE),
    REML2 = lmer(y ~ x + (x|g), weights = w, REML = TRUE),
    ML3 = lmer(y ~ x + z + (x+z|g), weights = w, REML = FALSE),
    REML3 = lmer(y ~ x + z + (x+z|g), weights = w, REML = TRUE))

compFunc <- function(lmeMod, lmerMod, tol = 1e-2){
    lmeOut <- c(as.numeric(VarCorr(lmeMod)[,"StdDev"]),
                as.numeric(summary(lmeMod)$tTable[,-c(3,5)]))

    lmerOut <- c(attr(VarCorr(lmerMod)$g, "stddev"),
                 attr(VarCorr(lmerMod), "sc"),
                 as.numeric(summary(lmerMod)$coefficients))
    names(lmeOut) <- names(lmerOut) <- NULL

    return(list(target = lmeOut, current = lmerOut, tolerance = tol))
}
comp <- mapply(compFunc, lmeMods, lmerMods, SIMPLIFY=FALSE)
stopifnot(all(sapply(comp, do.call, what = all.equal)))



## library("lme4")

## set.seed(2)
## n <- 40
## w <- runif(n)
## x <- runif(n)
## g <- factor(sample(1:10,n,replace=TRUE))
## Z <- model.matrix(~g-1);
## y <- Z%*%rnorm(ncol(Z)) + x + rnorm(n)/w^.5
## m <- lmer(y~x+(1|g),weights=w, REML = TRUE)

## fixef_lme4.0 <- c(-0.730654, 2.028954)
## stopifnot(all.equal(unname(fixef(m)), fixef_lme4.0, tol = 10^-3))

## sigma_lme4.0 <- 1.736143
## stopifnot(all.equal(sigma(m), sigma_lme4.0, tol = 10^-3))

## Sigma_lme4.0 <- 2.356705
## stopifnot(all.equal(as.vector(VarCorr(m)$g), Sigma_lme4.0, tol = 10^-3))

## SE_lme4.0 <- c(0.9507008, 1.3765086)
## stopifnot(all.equal(as.vector(summary(m)$coefficients[,2]), SE_lme4.0, tol = 10^-3))
