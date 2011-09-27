## Class definitions for the package

setClass("lmList",
         representation(call = "call",
                        pool = "logical"),
         contains = "list")

## TODO: export
setClass("lmList.confint", contains = "array")

### FIXME
## shouldn't we have "merPred"  with two *sub* classes "merPredD" and "merPredS"
## for the dense and sparse X cases ?

merPredD <-
setRefClass("merPredD", # Predictor class for mixed-effects models with dense X
            fields =
            list(
                 X = "ddenseModelMatrix",
                 Zt = "dgCMatrix",
                 Lambdat = "dgCMatrix",
                 Lind = "integer",
                 ptr = "externalptr",
                 beta0 = function(value)
                 if (missing(value)) .Call(merPredDbeta0, ptr) else .Call(merPredDsetBeta0, ptr, value),
                 theta = function(value)
                 if (missing(value)) .Call(merPredDtheta, ptr) else .Call(merPredDsetTheta, ptr, value),
                 u0 = function(value)
                 if (missing(value)) .Call(merPredDu0, ptr) else .Call(merPredDsetU0, ptr, value)
                 ),
            methods =
            list(
                 initialize = function(X, Zt, Lambdat, Lind, theta, ...) {
                     X <<- as(X, "ddenseModelMatrix")
                     Zt <<- as(Zt, "dgCMatrix")
                     Lambdat <<- as(Lambdat, "dgCMatrix")
                     Lind <<- as.integer(Lind)
                     stopifnot(all(sort(unique(Lind)) == seq_along(theta)))
                     ptr <<- .Call(merPredDCreate, X, Zt, Lambdat, Lind, theta)
                     callSuper(...)
                 },

                 CcNumer = function() {
                     'returns the numerator of the orthogonality convergence criterion'
                     .Call(merPredDCcNumer, ptr)
                 },
                 L = function() {
                     'returns the current value of the sparse Cholesky factor'
                     .Call(merPredDL, ptr)
                 },
                 LamtUt = function() {
                     'returns the current value of the product Lambdat %*% Ut'
                     .Call(merPredDLamtUt, ptr)
                 },
                 P = function() {
                     'returns the permutation vector for the sparse Cholesky factor'
                     .Call(merPredDPvec, ptr)
                 },
                 RX = function() {
                     'returns the dense downdated Cholesky factor for the fixed-effects parameters'
                     .Call(merPredDRX, ptr)
                 },
                 RXdiag = function() {
                     'returns the diagonal of the dense downdated Cholesky factor'
                     .Call(merPredDRXdiag, ptr)
                 },
                 RZX = function() {
                     'returns the cross term in Cholesky decomposition for all coefficients'
                     .Call(merPredDRZX, ptr)
                 },
                 Ut = function() {
                     'returns the transposed weighted random-effects model matrix'
                     .Call(merPredUt, ptr)
                 },
                 Utr = function() {
                     'returns the cross-product of the weighted random-effects model matrix\nand the weighted residuals'
                     .Call(merPredUtr, ptr)
                 },
                 V = function() {
                     'returns the weighted fixed-effects model matrix'
                     .Call(merPredDV, ptr)
                 },
                 VtV = function() {
                     'returns the weighted cross-product of the fixed-effects model matrix'
                     .Call(merPredDVtV, ptr)
                 },
                 Vtr = function() {
                     'returns the weighted cross-product of the fixed-effects model matrix\nand the residuals'
                     .Call(merPredVtr, ptr)
                 },
                 b = function(fac) {
                     'random effects on original scale for step factor fac'
                     .Call(merPredDb, ptr, as.numeric(fac))
                 },
                 beta = function(fac) {
                     'fixed-effects coefficients for step factor fac'
                     .Call(merPredDbeta, ptr, as.numeric(fac))
                 },
                 delb = function() {
                     'increment for the fixed-effects coefficients'
                     .Call(merPredDdelb, ptr)
                 },
                 delu = function() {
                     'increment for the orthogonal random-effects coefficients'
                     .Call(merPredDdelu, ptr)
                 },
                 getLambdat = function() {
                     'returns the current value of the relative covariance factor'
                     .Call(merPredDLambdat, ptr)
                 },
                 ldL2 = function() {
                     'twice the log determinant of the sparse Cholesky factor'
                     .Call(merPredDldL2, ptr)
                 },
                 ldRX2 = function() {
                     'twice the log determinant of the downdated dense Cholesky factor'
                     .Call(merPredDldRX2, ptr)
                 },
                 unsc = function() {
                     'the unscaled variance-covariance matrix of the fixed-effects parameters'
                     .Call(merPredDunsc, ptr)
                 },

                 linPred = function(fac) {
                     'evaluate the linear predictor for step factor fac'
                     .Call(merPredDlinPred, ptr, as.numeric(fac))
                 },
                 installPars = function(fac) {
                     'update u0 and beta0 to the values for step factor fac'
                     .Call(merPredDinstallPars, ptr, as.numeric(fac))
                 },
                 solve = function() {
                     'solve for the coefficient increments delu and delb'
                     .Call(merPredDsolve, ptr)
                 },
                 solveU = function() {
                     'solve for the coefficient increment delu only (beta is fixed)'
                     .Call(merPredDsolveU, ptr)
                 },
                 sqrL = function(fac) {
                     'squared length of u0 + fac * delu'
                     .Call(merPredDsqrL, ptr, as.numeric(fac))
                 },
                 u = function(fac) {
                     'orthogonal random effects for step factor fac'
                     .Call(merPredDu, ptr, as.numeric(fac))
                 },
                 updateDecomp = function() {
                     'update L, RZX and RX from Ut, Vt and VtV'
                     .Call(merPredDupdateDecomp, ptr)
                 },
                 updateRes = function(wtres) {
                     'update Vtr and Utr using the vector of weighted residuals'
                     .Call(merPredDupdateRes, ptr, as.numeric(wtres))
                 },
                 updateXwts = function(wts) {
                     'update Ut and V from Zt and X using X weights'
                     .Call(merPredDupdateXwts, ptr, as.numeric(wts))
                 }
                 )
            )
merPredD$lock("X", "Zt", "Lind")

setOldClass("family")

## FIXME --- glmerResp  &  lmerResp  should both extend a common (virtual) class
glmerResp <-
    setRefClass("glmerResp",
                fields =
                list(
                     family = "family",
                     y = "numeric",
                     ptr = "externalptr",
                     n = function(value) # only used in deviance residuals for binomial
                     if (missing(value)) .Call(glm_n, ptr) else .Call(glm_setN, ptr, as.numeric(value)),
                     offset = function(value)
                     if (missing(value)) .Call(lm_offset, ptr) else .Call(lm_setOffset, ptr, as.numeric(value)),
                     weights = function(value)
                     if (missing(value)) .Call(lm_weights, ptr) else .Call(lm_setWeights, ptr, as.numeric(value))
                     ),
                methods =
                list(
                     initialize = function(fam, y, ...) {
                         if (missing(y)) stop("response vector y must be specified")
                         y <<- as.numeric(y)
                         if (missing(fam) || !inherits(fam, "family"))
                             stop("glm family must be specified")
                         family <<- fam
                         ptr <<- .Call(glm_Create, family, y)
                         callSuper(...)
                     },
                     allInfo = function() {
                         'return all the information available on the object'
                         data.frame(y=y, n=n, offset=offset, weights=weights,
                                    eta=eta(), mu=mu(), muEta=muEta(), variance=variance(),
                                    sqrtrwt=sqrtrwt(), wtres=wtres(), sqrtXwt=sqrtXwt(),
                                    sqrtWrkWt=sqrtWrkWt(), wrkResids=wrkResids(),
                                    wrkResp=wrkResp())
                     },
                     devResid = function() {
                         'returns the vector of deviance residuals'
                         .Call(glm_devResid, ptr)
                     },
                     eta = function() {
                         'returns the linear predictor vector before any offset is added'
                         .Call(glm_eta, ptr)
                     },
                     fam = function() {
                         'returns the name of the glm family'
                         .Call(glm_family, ptr)
                     },
                     fitted = function() {
                         'returns the value of the linear predictor'
                         .Call(lm_mu, ptr)
                     },
                     Laplace = function(ldL2, ldRX2, sqrL) {
                         'returns the profiled deviance or REML criterion'
                         .Call(glm_Laplace, ptr, ldL2, ldRX2, sqrL)
                     },
                     link = function() {
                         'returns the name of the glm link'
                         .Call(glm_link, ptr)
                     },
                     mu = function() {
                         'returns the current mean response'
                         .Call(lm_mu, ptr)
                     },
                     muEta = function() {
                         'returns the diagonal of the Jacobian matrix, d mu/d eta'
                         .Call(glm_muEta, ptr)
                     },
                     resDev = function() {
                         'returns the sum of the deviance residuals'
                         .Call(glm_resDev, ptr)
                     },
                     sqrtXwt = function() {
                         'returns the square root of the X weights'
                         .Call(lm_sqrtXwt, ptr)
                     },
                     sqrtrwt = function() {
                         'returns the square root of the residual weights'
                         .Call(lm_sqrtrwt, ptr)
                     },
                     sqrtWrkWt = function() {
                         'returns the square root of the working X weights'
                         .Call(glm_sqrtWrkWt, ptr)
                     },
                     updateMu = function(gamma) {
                         'from the linear predictor, gamma, update the\nconditional mean response, mu, residuals and weights'
                         .Call(glm_updateMu, ptr, as.numeric(gamma))
                     },
                     updateWts = function() {
                         'update the residual and X weights from the current value of eta'
                         .Call(glm_updateWts, ptr)
                     },
                     variance = function() {
                         'returns the vector of variances'
                         .Call(glm_variance, ptr)
                     },
                     wrkResids = function() {
                         'returns the vector of working residuals'
                         .Call(glm_wrkResids, ptr)
                     },
                     wrkResp = function() {
                         'returns the vector of working residuals'
                         .Call(glm_wrkResp, ptr)
                     },
                     wrss = function() {
                         'returns the weighted residual sum of squares'
                         .Call(lm_wrss, ptr)
                     },
                     wtres = function() {
                         'returns the vector of weighted residuals'
                         .Call(lm_wtres, ptr)
                     })
                )

glmerResp$lock("family", "y")

## seems currently *unused* -FIXME-
glmFamily <-
    setRefClass("glmFamily",
                fields =
                list(
                     family = "family",
                     ptr    = "externalptr"
                     ),
                methods =
                list(
                     initialize = function(fam, ...) {
                         stopifnot(inherits(fam, "family"))
                         family <<- fam
                         ptr <<- .Call(glmFamily_Create, fam)
                         callSuper(...)
                     },
                     link = function(mu) {
                         'applies the (forward) link function to mu'
                         .Call(glmFamily_link, ptr, as.numeric(mu))
                     },
                     linkInv = function(eta) {
                         'applies the inverse link function to eta'
                         .Call(glmFamily_linkInv, ptr, as.numeric(eta))
                     },
                     muEta = function(eta) {
                         'applies the muEta function to eta'
                         .Call(glmFamily_muEta, ptr, as.numeric(eta))
                     },
                     variance = function(mu) {
                         'applies the variance function to mu'
                         .Call(glmFamily_variance, ptr, as.numeric(mu))
                     },
                     devResid = function(mu, weights, y) {
                         'applies the devResid function to mu, weights and y'
                         mu <- as.numeric(mu)
                         weights <- as.numeric(weights)
                         y <- as.numeric(y)
                         stopifnot(length(mu) == length(weights),
                                   length(mu) == length(y),
                                   all(weights >= 0))
                         .Call(glmFamily_devResid, ptr, mu, weights, y)
                     }
                 )
                )

lmerResp <-
    setRefClass("lmerResp",
                fields =
                list(
                     y = "numeric",
                     ptr = "externalptr",
                     offset = function(value)
                     if (missing(value)) .Call(lm_offset, ptr) else .Call(lm_setOffset, ptr, as.numeric(value)),
                     weights = function(value)
                     if (missing(value)) .Call(lm_weights, ptr) else .Call(lm_setWeights, ptr, as.numeric(value)),
                     REML = function(value)
                     if (missing(value)) .Call(lmer_REML, ptr) else .Call(lmer_setREML, ptr, value)
                     ),
                methods =
                list(
                     initialize = function(y, ...) {
                         if (missing(y)) stop("response vector y must be specified")
                         y <<- as.numeric(y)
                         ptr <<- .Call(lmer_Create, y)
                         callSuper(...)
                     },
                     fitted = function() {
                         'returns the value of the linear predictor'
                         .Call(lm_mu, ptr)
                     },
                     Laplace = function(ldL2, ldRX2, sqrL) {
                         'returns the profiled deviance or REML criterion'
                         .Call(lmer_Laplace, ptr, ldL2, ldRX2, sqrL)
                     },
                     updateMu = function(gamma) {
                         'from the linear predictor, gamma, update the\nconditional mean response, mu, residuals and weights'
                         .Call(lm_updateMu, ptr, as.numeric(gamma))
                     },
                     wrss = function() {
                         'returns the weighted residual sum of squares'
                         .Call(lm_wrss, ptr)
                     },
                     wtres = function() {
                         'returns the vector of weighted residuals'
                         .Call(lm_wtres, ptr)
                     })
                )

lmerResp$lock("y")

setClass("merMod",
         representation(call    = "call",
			frame   = "data.frame", # "model.frame" is not S4-ized yet
                        flist   = "list",
                        cnms    = "list",
                        lower   = "numeric",
                        Gp      = "integer",
                        theta   = "numeric",
                        beta    = "numeric",
                        u       = "numeric",
                        devcomp = "list",
                        pp      = "merPredD"))

setClass("lmerMod", representation(resp = "lmerResp"), contains = "merMod")

setClass("glmerMod", representation(resp = "glmerResp"), contains = "merMod")
