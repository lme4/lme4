######################################################################
## Pure R implementation of GLMM deviances (not meant to be efficient,
## but rather for sanity checks)
##
## Based on the following set of distinctions:
##
##            conditional          unconditional
## relative   GLM deviance         ???
## absolute   -2 GLM log density   standard normal average of -2 GLM density
##
## Naming conventions:
##
## rel vs abs: 'rel'ative to the deviance of the saturated model, or
## 'abs'olute?
##
## cond vs uncond: 'cond'itional on values for the spherical random
## effects, or 'uncond'itional?
##
## residuals*: conditional deviance residuals (unconditional deviance
## residuals do not really make sense to me), which are either
## relative or absolute.
##
## devianceCond*: sum(residuals*)
##
## devianceUncondRel: does not make sense to me, so just throws an
## error.
##
## devianceUncondAbs: requires a difficult integral, and so calls one
## of the following three functions that involve different numerical
## approximations to this integral: (1) devianceLaplace, (2)
## devianceAGQ, or (3) devianceImportSamp
##
######################################################################

library(lme4)
library(mvtnorm)

######################################################################
##' Two types of conditional deviance residuals
##' 
##' @param object merMod object
##' @param lp optional vector of linear predictors
residualsRel <- function(object, lp) {
    ## family()$dev.resid
    if(missing(lp)) return(object@resp$devResid())
    fam <- family(object)
    resFunc <- fam$dev.resid
    resFunc(y = object@resp$y,
            mu = fam$linkinv(lp),
            wt = weights(object))
}
residualsAbs <- function(object, lp, lpForScale) {
    ## family()$aic but with three modifications:
    ## (1) use mapply to compute separately for each observation.
    ## (2) the devFamAIC wrapper is a hack to deal with estimation of
    ## scale parameters (when they exist).
    ## (3) -2*useSc term removes the AIC complexity penalty on the
    ## scale parameter (which is present in family objects even though
    ## complexity from regression coefficients is completly absent).
    if(missing(lp)) lp <- condLinPred(object)
    fam <- family(object)
    aicFunc <- fam$aic
    useSc <- getME(object, "devcomp")$dims["useSc"]
    if((useSc == 1) && (missing(lpForScale))) lpForScale <- lp
    mapply(aicFunc,
           y = object@resp$y,
           n = object@resp$n,
           mu = fam$linkinv(lp),
           wt = weights(object),
           dev = devFamAIC(object, lpForScale)) - 2 * useSc
}

######################################################################
##' Four types of deviance
##'
##' @param object merMod object
##' @param lp optional linear predictor
##' @param approximation type of integral approximation
##' @param NAGQ number of quadrature points
##' @param Nsim number of importance sampling draws
devianceCondRel <- function(object, lp            ) sum(residualsRel(object, lp            ))
devianceCondAbs <- function(object, lp, lpForScale) sum(residualsAbs(object, lp, lpForScale))
devianceUncondRel <- function(object) stop("relative-unconditional deviance not defined")
devianceUncondAbs <- function(object,
                              approximation = c("Laplace", "AGQ", "importSamp"),
                              NAGQ = 5, Nsim = 100) {
    switch(approximation,
           Laplace    = devianceLaplace   (object),
           AGQ        = devianceAGQ       (object, NAGQ),
           importSamp = devianceImportSamp(object, Nsim))
}

######################################################################
##' Three approaches to approximating the unconditional absolute
##' deviance: devianceLaplace, devianceAGQ, and devianceImportSamp
devianceLaplace <- function(object) {
    ldL2(object) + sqrL(object) + devianceCondAbs(object)
}
devianceAGQ <- function(object, NAGQ) {
    if(length(facs <- object@flist) != 1L) stop("only single scalar random effects allowed for AGQ")
    fac <- facs[[1]]
    q <- getME(object, "q")
    AGQinfo <- setupAGQ(object, NAGQ)
    lp <- apply(AGQinfo$u, 2, condLinPred, object = object)
    resAbs <- apply(lp, 2, residualsAbs, object = object, lpForScale = condLinPred(object))
    logCondDenQuad <- -0.5 * apply(resAbs, 2, tapply, fac, "sum")
    ## FIXME: unstable for large sample sizes (logCondDenQuad's too
    ## large in those cases), so don't do a matrix product but
    ## something that keeps things on the log scale for longer
    summandAGQ <- c(exp(logCondDenQuad - 0.5*(AGQinfo$u^2)) %*% with(AGQinfo, w * exp(-l)))
    ldL2(object) + q * log(2 * pi) - 2 * sum(log(summandAGQ))
}
devianceImportSamp <- function(object, Nsim) {
    require(mvtnorm)
    Iq <- Diagonal(getME(object, "q"))
                                        # conditional covariance
                                        # matrix for spherical random
                                        # effects
    uCondVar <- as.matrix(crossprod(Iq, solve(getME(object, "L"), Iq, system = "A")))
                                        # simulated spherical random
                                        # effects
    uSim <- t(rmvnorm(Nsim, getME(object, "u"), sigma = uCondVar))
    lp <- apply(uSim, 2, condLinPred, object = object)
        condDen <- -0.5 * apply(lp, 2, devianceCondAbs, object =
                                object, lpForScale = condLinPred(object))
                                        # target density of spherical
                                        # random effects
      targetDen <- apply(dnorm(uSim, log = TRUE), 2, sum)
                                        # importance sampling density
    samplingDen <- apply(uSim, 2, dmvnorm, getME(object, "u"), uCondVar, TRUE)
                                        # importance sampling estimate
    logIntegrands <- condDen + targetDen - samplingDen
    -2 * (log(mean(exp(logIntegrands - max(logIntegrands)))) + max(logIntegrands))
}

######################################################################
##' helpers
##'
##' @param object merMod object
##' @param u q-vector of spherical random effects
condLinPred <- function(object, u) {
    ## conditional linear predictor
    if(missing(u)) u <- getME(object, "u")
    fe <- getME(object, "X") %*% getME(object, "beta")
    re <- getME(object, "Z") %*% getME(object, "Lambda") %*% u
    off <- getME(object, "offset")
    as.numeric(fe + re + off)
}
setupAGQ <- function(object, NAGQ) {
    ## setup an AGQ grid
    if(missing(NAGQ)) NAGQ <- getME(object, "devcomp")$dims["nAGQ"]
    if(NAGQ < 2) stop("must have at least two quadrature nodes (NAGQ > 1)")
    ghr <- GHrule(NAGQ)
    z <- ghr[,"z"] # quadrature nodes
    w <- ghr[,"w"] # quadrature weights
    l <- ghr[,"ldnorm"] # dnorm(z, log = TRUE)
    sig <- 1/getME(object, "L")@x # scale for the change of variables
    list(u = sweep(outer(sig, z), 1, getME(object, "u"), "+"),
         z = z, w = w, l = l)
}
sqrL <- function(object) {
    ## squared length
    if(isGLMM(object)) return(object@pp$sqrL(1))
    if(is.numeric(object)) return(sum(object^2))
    stop("cannot compute squared length for this kind of object")
}
ldL2 <- function(object) {
    ## twice the log-determinant of the sparse Cholesky factor
    object@pp$ldL2()
}
devFamAIC <- function(object, lp) {
    ## hack to deal with scale parameters
    fam <- family(object)$family
    n <- getME(object, "n")
    devPerObs <- rep(devianceCondRel(object, lp)/n, n)
    devPerWts <- weights(object) * (devianceCondRel(object, lp)/sum(weights(object)))
    switch(fam,
           gaussian         = devPerObs,
           Gamma            = devPerWts,
           inverse.gaussian = devPerWts,
           rep(NA, n)) ## default case without scale parameters
}


######################################################################
## examples
######################################################################

######################################################################
## basic one (looks good)
set.seed(1)
gm <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
            data = cbpp, family = binomial)
dqu <- c(devianceLaplace(gm), sapply(2:25, devianceAGQ, object = gm))
Nis <- 2000
dis <- devianceImportSamp(gm, Nis)
plot(1:25, dqu,
     type = "o", xlab = "number of quadrature points", ylab = "deviance",
     ylim = range(dqu, dis))
abline(h = dis, col = "blue", lty = 2) # will be some randomness so
                                       # don't expect perfect
                                       # congruence with asymptotic
                                       # AGQ results (increase Nis for
                                       # more accuracy (hopefully!))
all.equal(deviance(gm), devianceCondRel(gm))
all.equal(-2 * c(logLik(gm)), devianceLaplace(gm))

######################################################################
## grouseticks (looks good)
set.seed(1)
data(grouseticks)
form <- TICKS ~ YEAR + HEIGHT + (1 | BROOD) ## + (1 | INDEX) + (1 | LOCATION)
gm  <- glmer(form, family = "poisson", data = grouseticks)
dqu <- c(devianceLaplace(gm), sapply(2:25, devianceAGQ, object = gm))
dis <- devianceImportSamp(gm, Nis)
plot(1:25, dqu,
     type = "o", xlab = "number of quadrature points", ylab = "deviance",
     ylim = range(dqu, dis))
abline(h = dis, col = "blue", lty = 2)
all.equal(deviance(gm), devianceCondRel(gm))
all.equal(-2 * c(logLik(gm)), devianceLaplace(gm))

######################################################################
## based on tests/glmmExt.R with Gamma (looks good)
set.seed(1)
d <- expand.grid(block=LETTERS[1:26], rep=1:10, KEEP.OUT.ATTRS = FALSE)
d$x <- runif(nrow(d))  ## sd=1
reff_f <- rnorm(length(levels(d$block)),sd=1)
## need intercept large enough to avoid negative values
d$eta0 <- 4 + 3 * d$x  ## fixed effects only
d$eta <- d$eta0 + reff_f[d$block]
dgl <- d
dgl$mu <- exp(d$eta)
dgl$y <- rgamma(nrow(d), scale = dgl$mu/2, shape = 2)
gm <- glmer(y ~ 1 + (1 | block), data = dgl, family = Gamma(link = "log"))
dqu <- c(devianceLaplace(gm), sapply(2:25, devianceAGQ, object = gm))
dis <- devianceImportSamp(gm, Nis)
plot(1:25, dqu,
     type = "o", xlab = "number of quadrature points", ylab = "deviance",
     ylim = range(dqu, dis))
abline(h = dis, col = "blue", lty = 2)
all.equal(deviance(gm), devianceCondRel(gm))
all.equal(-2 * c(logLik(gm)) - 2, ## minus two is for the scale
                                  ## parameter (should include this in
                                  ## lme4)
          devianceLaplace(gm))

aic1 <- Gamma()$aic(getME(gm, "y"),
                    gm@resp$n,
                    getME(gm, "mu"),
                    weights(gm),
                    deviance(gm))

aic0 <- Gamma()$aic(getME(gm, "y"),
                    gm@resp$n,
                    getME(gm, "y"),
                    weights(gm),
                    deviance(gm))

disp <- deviance(gm)/sum(weights(gm))

aic1 - aic0
sum(Gamma()$dev.resids(getME(gm, "y"),
                       getME(gm, "mu"),
                       weights(gm)))/disp



