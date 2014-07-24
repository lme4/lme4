compFunc <- function(lmeMod, lmerMod, tol = 1e-2){
    lmeVarCorr <- nlme:::VarCorr(lmeMod)[,"StdDev"]
    lmeCoef <- summary(lmeMod)$tTable[,-c(3,5)]
    lmeOut <- c(as.numeric(lmeVarCorr), as.numeric(lmeCoef))
    keep <- !is.na(lmeOut)
    lmeOut <- lmeOut[keep]
    dn <- dimnames(lmeCoef)
    if(is.null(dn)) dn <- list("", names(lmeCoef))
    names(lmeOut) <- c(
        paste(names(lmeVarCorr), "Var"),
        as.character(do.call(outer, c(dn, list("paste")))))[keep]

    ## get nested RE variances in the same order as nlme
    ## FIXME:  not sure if this works generally
    vcLmer <- VarCorr(lmerMod)
    vcLmer <- vcLmer[length(vcLmer):1]
    ##

    lmerVarCorr <- c(sapply(vcLmer, attr, "stddev"),
                     attr(VarCorr(lmerMod), "sc"))
    ## differentiate lme4{new} and lme4.0 :
    lmerCoef <- if(is(lmerMod, "merMod"))
	summary(lmerMod)$coefficients else summary(lmerMod)@coefs
    lmerOut <- c(lmerVarCorr, as.numeric(lmerCoef))
    names(lmerOut) <- names(lmeOut)

    return(list(target = lmeOut, current = lmerOut, tolerance = tol))
}

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
residSD <- rpois(nObs,1) + 1
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
    ML1   = lme(y ~ x,   random = ~ 1|g,   weights = varFixed(~v), method = "ML"),
    REML1 = lme(y ~ x,   random = ~ 1|g,   weights = varFixed(~v), method = "REML"),
    ML2   = lme(y ~ x,   random = ~ x|g,   weights = varFixed(~v), method = "ML"),
    REML2 = lme(y ~ x,   random = ~ x|g,   weights = varFixed(~v), method = "REML"),
    ML1   = lme(y ~ x+z, random = ~ x+z|g, weights = varFixed(~v), method = "ML"),
    REML2 = lme(y ~ x+z, random = ~ x+z|g, weights = varFixed(~v), method = "REML"))

library("lme4")
lmerMods <- list(
    ML1 =   lmer(y ~ x +  (1|g),      weights = w, REML = FALSE),
    REML1 = lmer(y ~ x +  (1|g),      weights = w, REML = TRUE),
    ML2 =   lmer(y ~ x +  (x|g),      weights = w, REML = FALSE),
    REML2 = lmer(y ~ x +  (x|g),      weights = w, REML = TRUE),
    ML3 =   lmer(y ~ x + z + (x+z|g), weights = w, REML = FALSE),
    REML3 = lmer(y ~ x + z + (x+z|g), weights = w, REML = TRUE))

comp <- mapply(compFunc, lmeMods, lmerMods, SIMPLIFY=FALSE)
stopifnot(all(sapply(comp, do.call, what = all.equal)))
## Look at the relative differences:
sapply(mapply(compFunc, lmeMods, lmerMods, SIMPLIFY=FALSE, tol = 0),
       do.call, what = all.equal)

## add simulated weights to the sleepstudy example
n <- nrow(sleepstudy)
v <- rpois(n,1) + 1
w <- 1/v
sleepLme <- lme(Reaction ~ Days, random = ~ Days|Subject,
                sleepstudy, weights = varFixed(~v),
                method = "ML")
sleepLmer <- lmer(Reaction ~ Days + (Days|Subject),
                  sleepstudy, weights = w,
                  REML = FALSE)
sleepComp <- compFunc(sleepLme, sleepLmer)
stopifnot(do.call(all.equal, sleepComp))
## look at relative differences:
sleepComp$tolerance <- 0
do.call(all.equal, sleepComp)

if (require("mlmRev")) {
    n <- nrow(Chem97)
    v <- rpois(n,1) + 1
    w <- 1/v
    Chem97Lme <- lme(score ~ 1, random = ~ 1|lea/school, Chem97)
    Chem97Lmer <- lmer(score ~ (1|lea/school), Chem97)
    Chem97Comp <- compFunc(Chem97Lme, Chem97Lmer)
    stopifnot(do.call(all.equal, Chem97Comp))
    ## look at relative differences:
    Chem97Comp$tolerance <- 0
    do.call(all.equal, Chem97Comp)
}

set.seed(2)
n <- 40
w <- runif(n)
x <- runif(n)
g <- factor(sample(1:10,n,replace=TRUE))
Z <- model.matrix(~g-1);
y <- Z%*%rnorm(ncol(Z)) + x + rnorm(n)/w^.5
m <- lmer(y ~ x + (1|g), weights=w, REML = TRUE)

## CRAN-forbidden:
## has4.0 <- require("lme4.0"))
has4.0 <- FALSE
if(has4.0) {
    ## m.0 <- lme4.0::lmer(y ~ x + (1|g), weights=w, REML = TRUE)
    lmer0 <- get("lmer", envir=asNamespace("lme4.0"))
    m.0 <- lmer0(y ~ x + (1|g), weights=w, REML = TRUE)
    dput(fixef(m.0)) # c(-0.73065400610675, 2.02895402562926)
    dput(sigma(m.0)) # 1.73614301673377
    dput(VarCorr(m.0)$g[1,1]) # 2.35670451590395
    dput(unname(coef(summary(m.0))[,"Std. Error"]))
    ## c(0.95070076853232, 1.37650858268602)
}
fixef_lme4.0 <- c(-0.7306540061, 2.0289540256)
sigma_lme4.0 <- 1.7361430
Sigma_lme4.0 <- 2.3567045
SE_lme4.0 <- c(0.95070077, 1.37650858)
if(has4.0) try(detach("package:lme4.0"))

stopifnot(all.equal(unname(fixef(m)), fixef_lme4.0, tolerance = 1e-3))
          all.equal(unname(fixef(m)), fixef_lme4.0, tolerance = 0) #-> 1.657e-5

## but these are not at all equal :
(all.equal(sigma(m),                		sigma_lme4.0, tolerance = 10^-3)) # 0.4276
(all.equal(as.vector(VarCorr(m)$g), 		Sigma_lme4.0, tolerance = 10^-3)) # 1.038
(all.equal(as.vector(summary(m)$coefficients[,2]), SE_lme4.0, tolerance = 10^-3)) # 0.4276
## so, lme4.0 was clearly wrong here


##' make sure models that differ only in a constant
##' prior weight have identical deviance:
fm <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy,REML=FALSE)
fm_wt <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, weights = rep(5, nrow(sleepstudy)),REML=FALSE)
all.equal(deviance(fm), deviance(fm_wt))
