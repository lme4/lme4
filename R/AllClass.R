## Class definitions for the package

setClass("lmList",
         representation(call = "call",
                        pool = "logical"),
         contains = "list")

## TODO: export
setClass("lmList.confint", contains = "array")

merPredD <-                 # Do we need the generator object? Probably - at least while debugging
    setRefClass("merPredD", # Predictor class for mixed-effects models with dense X
                fields =
                list(
                     X = "ddenseModelMatrix",
                     Z = "dgCMatrix",
                     Lambda = "dgCMatrix",
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
                     initialize = function(...) {
                         dots <- list(...)
                         stopifnot(is(dots$X, "ddenseModelMatrix"),
                                   is(dots$Z, "dgCMatrix"),
                                   is(dots$Lambda, "CsparseMatrix"),
                                   !is.null(dots$Lind),
                                   !is.null(dots$theta))
                         X <<- dots$X
                         Z <<- dots$Z
                         Lambda <<- dots$Lambda
                         Lind <<- as.integer(dots$Lind)
                         stopifnot(all(sort(unique(Lind)) == seq_along(dots$theta)))
                         ptr <<- .Call(merPredDCreate, X, Z, Lambda, Lind, dots$theta)
                         if (!is.null(dots$beta0)) .Call(merPredDsetBeta0, ptr, dots$beta0)
                         if (!is.null(dots$u0)) .Call(merPredDsetBeta0, ptr, dots$u0)
                     },
                     getLambda = function() .Call(merPredDLambda, ptr),
                     RX = function() .Call(merPredDRX, ptr),
                     RZX = function() .Call(merPredDRZX, ptr),
                     delb = function() .Call(merPredDdelb, ptr),
                     delu = function() .Call(merPredDdelu, ptr),
                     ldL2 = function() .Call(merPredDldL2, ptr),
                     ldRX2 = function() .Call(merPredDldRX2, ptr),
                     getPvec = function() .Call(merPredDPvec, ptr),
                     unsc = function() .Call(merPredDunsc, ptr),

                     linPred = function(fac) .Call(merPredDlinPred, ptr, as.numeric(fac)),
                     installPars = function(fac) .Call(merPredDinstallPars, ptr, as.numeric(fac)),
                     sqrL = function(fac) .Call(merPredDsqrL, ptr, as.numeric(fac))
                     )
                )
merPredD$lock("X", "Z", "Lind")

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
                     initialize = function(...) {
                         dots <- list(...)
                         if (is.null(dots$y)) stop("response vector y must be specified")
                         y <<- as.numeric(dots$y)
                         ptr <<- .Call(lmerRespCreate, y)
                         if (!is.null(dots$REML)) REML(as.integer(dots$REML)[1])
                         if (length(dots$offset)) offset(as.numeric(dots$offset))
                         if (length(dots$weights)) weights(as.numeric(dots$weights))
                     },
                     wrss = function() {
                         'returns the weighted residual sum of squares'
                         .Call(modRespwrss, ptr)
                     },
                     wtres = function() {
                         'returns the vector of weighted residuals'
                         .Call(modRespwtres, ptr)
                     },
                     Laplace = function(ldL2, ldRX2, sqrL) {
                         'returns the profiled deviance or REML criterion'
                         .Call(lmerRespLaplace, ptr, ldL2, ldRX2, sqrL)
                     },
                     updateMu = function(gamma) {
                         'from the linear predictor, gamma, update the\nconditional mean response, mu, residuals and weights'
                         .Call(lmerRespupdateMu, ptr, as.numeric(gamma))
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
