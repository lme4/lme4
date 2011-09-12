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

merPredD <-                 # Do we need the generator object? Probably - at least while debugging
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
                 L = function() {
                     'returns the current value of the sparse Cholesky factor'
                     .Call(merPredDL, ptr)
                 },
                 P = function() {
                     'returns the permutation vector for the sparse Cholesky factor'
                     .Call(merPredDPvec, ptr)
                 },
                 RX = function() {
                     'returns the dense downdated Cholesky factor for the fixed-effects parameters'
                     .Call(merPredDRX, ptr)
                 },
                 RXdiag = function() .Call(merPredDRXdiag, ptr),
                 RZX = function() .Call(merPredDRZX, ptr),
                 VtV = function() .Call(merPredDVtV, ptr),
                 delb = function() .Call(merPredDdelb, ptr),
                 delu = function() .Call(merPredDdelu, ptr),
                 getLambda = function() .Call(merPredDLambda, ptr),
                 ldL2 = function() .Call(merPredDldL2, ptr),
                 ldRX2 = function() .Call(merPredDldRX2, ptr),
                 unsc = function() .Call(merPredDunsc, ptr),

                 linPred = function(fac) .Call(merPredDlinPred, ptr, as.numeric(fac)),
                 installPars = function(fac) .Call(merPredDinstallPars, ptr, as.numeric(fac)),
                 sqrL = function(fac) .Call(merPredDsqrL, ptr, as.numeric(fac))
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
                     offset = function(value)
                     if (missing(value)) .Call(modRespoffset, ptr) else .Call(modRespsetOffset, ptr, as.numeric(value)),
                     weights = function(value)
                     if (missing(value)) .Call(modRespweights, ptr) else .Call(modRespsetWeights, ptr, as.numeric(value))
                     ),
                methods =
                list(
                     initialize = function(fam, y, ...) {
                         if (missing(y)) stop("response vector y must be specified")
                         y <<- as.numeric(y)
                         if (missing(fam) || !inherits(fam, "family"))
                             stop("glm family must be specified")
                         family <<- fam
                         ptr <<- .Call(glmerRespCreate, family, y)
                         callSuper(...)
                     },
                     devResid = function() {
                         'returns the vector of deviance residuals'
                         .Call(glmerRespdevResid, ptr)
                     },
                     eta = function() {
                         'returns the linear predictor vector before any offset is added'
                         .Call(glmerRespeta, ptr)
                     },
                     fam = function() {
                         'returns the name of the glm family'
                         .Call(glmerRespfamily, ptr)
                     },
                     fitted = function() {
                         'returns the value of the linear predictor'
                         .Call(modRespmu, ptr)
                     },
                     Laplace = function(ldL2, ldRX2, sqrL) {
                         'returns the profiled deviance or REML criterion'
                         .Call(glmerRespLaplace, ptr, ldL2, ldRX2, sqrL)
                     },
                     link = function() {
                         'returns the name of the glm link'
                         .Call(glmerResplink, ptr)
                     },
                     updateMu = function(gamma) {
                         'from the linear predictor, gamma, update the\nconditional mean response, mu, residuals and weights'
                         .Call(glmerRespupdateMu, ptr, as.numeric(gamma))
                     },
                     wrkResid = function() {
                         'returns the vector of working residuals'
                         .Call(modRespwrkResid, ptr)
                     },
                     wrkResp = function() {
                         'returns the vector of working residuals'
                         .Call(modRespwrkResp, ptr)
                     },
                     wrss = function() {
                         'returns the weighted residual sum of squares'
                         .Call(modRespwrss, ptr)
                     },
                     wtres = function() {
                         'returns the vector of weighted residuals'
                         .Call(modRespwtres, ptr)
                     })
                )

glmerResp$lock("family", "y")

lmerResp <-
    setRefClass("lmerResp",
                fields =
                list(
                     y = "numeric",
                     ptr = "externalptr",
                     offset = function(value)
                     if (missing(value)) .Call(modRespoffset, ptr) else .Call(modRespsetOffset, ptr, as.numeric(value)),
                     weights = function(value)
                     if (missing(value)) .Call(modRespweights, ptr) else .Call(modRespsetWeights, ptr, as.numeric(value)),
                     REML = function(value)
                     if (missing(value)) .Call(lmerRespREML, ptr) else .Call(lmerRespsetREML, ptr, value)
                     ),
                methods =
                list(
                     initialize = function(y, ...) {
                         if (missing(y)) stop("response vector y must be specified")
                         y <<- as.numeric(y)
                         ptr <<- .Call(lmerRespCreate, y)
                         callSuper(...)
                     },
                     fitted = function() {
                         'returns the value of the linear predictor'
                         .Call(modRespmu, ptr)
                     },
                     Laplace = function(ldL2, ldRX2, sqrL) {
                         'returns the profiled deviance or REML criterion'
                         .Call(lmerRespLaplace, ptr, ldL2, ldRX2, sqrL)
                     },
                     updateMu = function(gamma) {
                         'from the linear predictor, gamma, update the\nconditional mean response, mu, residuals and weights'
                         .Call(lmerRespupdateMu, ptr, as.numeric(gamma))
                     },
                     wrss = function() {
                         'returns the weighted residual sum of squares'
                         .Call(modRespwrss, ptr)
                     },
                     wtres = function() {
                         'returns the vector of weighted residuals'
                         .Call(modRespwtres, ptr)
                     })
                )

lmerResp$lock("y")


setClass("merMod",
         representation(call    = "call",
			frame   = "data.frame", # "model.frame" is not S4-ized yet
                        flist   = "list",
                        cnms    = "list",
                        lower   = "numeric",
                        theta   = "numeric",
                        beta    = "numeric",
                        u       = "numeric",
                        devcomp = "list",
                        pp      = "merPredD"))

setClass("lmerMod", representation(resp = "lmerResp"), contains = "merMod")
