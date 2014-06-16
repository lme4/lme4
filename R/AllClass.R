### Class definitions for the package

##' Class "lmList" of 'lm' Objects on Common Model
##'
##' Class \code{"lmList"} is an S4 class with basically a list of objects of
##' class \code{\link{lm}} with a common model.
##' @name lmList-class
##' @aliases lmList-class show,lmList-method
##' @docType class
##' @section Objects from the Class: Objects can be created by calls of the form
##' \code{new("lmList", ...)} or, more commonly, by a call to
##' \code{\link{lmList}}.
##' @keywords classes
##' @export
setClass("lmList",
         representation(call = "call",
                        pool = "logical"),
         contains = "list")

## TODO: export
setClass("lmList.confint", contains = "array")

forceCopy <- function(x) {
    if (is.numeric(x)) return(x+0)
    ## FIXME: doesn't handle non-numeric fields yet ...
    ## resp@family is the only non-numeric field other than Ptr
    ## need equivalent force-copy no-op for other classes
    x
}

### FIXME
### shouldn't we have "merPred"  with two *sub* classes "merPredD" and "merPredS"
### for the dense and sparse X cases ?

##' Generator object for the \code{\linkS4class{merPredD}} class
##'
##' The generator object for the \code{\linkS4class{merPredD}} reference class.
##' Such an object is primarily used through its \code{new} method.
##'
##' @param ... List of arguments (see Note).
##' @section Methods:
##' \describe{
##'   \item{new(X, Zt, Lambdat, Lind, theta, n):}{Create a new \code{\linkS4class{merPredD}} object}
##' }
##' @note Arguments to the \code{new} methods must be named arguments:
##' \itemize{
##' \item{X}{ dense model matrix for the fixed-effects parameters, to be stored
##'     in the \code{X} field.}
##' \item{Zt}{ transpose of the sparse model matrix for the random effects.  It
##'     is stored in the \code{Zt} field.}
##' \item{Lambdat}{ transpose of the sparse lower triangular relative variance
##'     factor (stored in the \code{Lambdat} field).}
##' \item{Lind}{ integer vector of the same length as the \code{"x"} slot in the
##'     \code{Lambdat} field.  Its elements should be in the range 1 to the length
##'     of the \code{theta} field.}
##' \item{theta}{ numeric vector of variance component parameters (stored in the
##'     \code{theta} field).}
##' \item{n}{ sample size, usually \code{nrow(X)}.}
##' }
##' @seealso \code{\linkS4class{merPredD}}
##' @keywords classes
##' @export
merPredD <-
    setRefClass("merPredD", # Predictor class for mixed-effects models with dense X
                fields =
                list(Lambdat = "dgCMatrix",   # depends: theta and Lind
                     LamtUt  = "dgCMatrix",   # depends: Lambdat and Ut
                     Lind    = "integer",     # depends: nothing
                     Ptr     = "externalptr", # depends: 
                     RZX     = "matrix",      # depends: lots
                     Ut      = "dgCMatrix",   # depends: Zt and weights
                     Utr     = "numeric",     # depends: lots
                     V       = "matrix",      # depends: 
                     VtV     = "matrix",
                     Vtr     = "numeric",
                     X       = "matrix",
                     Xwts    = "numeric",
                     Zt      = "dgCMatrix",
                     beta0   = "numeric",
                     delb    = "numeric",
                     delu    = "numeric",
                     theta   = "numeric",
                     u0      = "numeric"),
                methods =
                list(
                     initialize = function(X, Zt, Lambdat, Lind, theta, n, ...) {
                         if (!nargs()) return
                         X <<- as(X, "matrix")
                         Zt <<- as(Zt, "dgCMatrix")
                         Lambdat <<- as(Lambdat, "dgCMatrix")
                         Lind <<- as.integer(Lind)
                         theta <<- as.numeric(theta)
                         N <- nrow(X)
                         p <- ncol(X)
                         q <- nrow(Zt)
                         stopifnot(length(theta) > 0L,
                                   length(Lind) > 0L,
                                   all(sort(unique(Lind)) == seq_along(theta)))
                         RZX <<- array(0, c(q, p))
                         Utr <<- numeric(q)
                         V <<- array(0, c(n, p))
                         VtV <<- array(0, c(p, p))
                         Vtr <<- numeric(p)
                         b0 <- list(...)$beta0
                         beta0 <<- if (is.null(b0)) numeric(p) else b0
                         delb <<- numeric(p)
                         delu <<- numeric(q)
                         uu <- list(...)$u0
                         u0 <<- if (is.null(uu)) numeric(q) else uu
                         Ut <<- if (n == N) Zt + 0 else
                             Zt %*% sparseMatrix(i=seq_len(N), j=as.integer(gl(n, 1, N)), x=rep.int(1,N))
                         ## The following is a kludge to overcome problems when Zt is square
                         ## by making LamtUt rectangular
                         LtUt <- Lambdat %*% Ut
                         ## if (nrow(LtUt) == ncol(LtUt))
                         ##     LtUt <- cbind2(LtUt,
                         ##                    sparseMatrix(i=integer(0),
                         ##                                 j=integer(0),
                         ##                                 x=numeric(0),
                         ##                                 dims=c(nrow(LtUt),1)))
                         LamtUt <<- LtUt
                         Xw <- list(...)$Xwts
                         Xwts <<- if (is.null(Xw)) rep.int(1, N) else as.numeric(Xw)
                         initializePtr()
                     },
                     CcNumer      = function() {
                         'returns the numerator of the orthogonality convergence criterion'
                         .Call(merPredDCcNumer, ptr())
                     },
                     L            = function() {
                         'returns the current value of the sparse Cholesky factor'
                         .Call(merPredDL, ptr())
                     },
                     P            = function() {
                         'returns the permutation vector for the sparse Cholesky factor'
                         .Call(merPredDPvec, ptr())
                     },
                     RX           = function() {
                         'returns the dense downdated Cholesky factor for the fixed-effects parameters'
                         .Call(merPredDRX, ptr())
                     },
                     RXi          = function() {
                         'returns the inverse of the dense downdated Cholesky factor for the fixed-effects parameters'
                         .Call(merPredDRXi, ptr())
                     },
                     RXdiag       = function() {
                         'returns the diagonal of the dense downdated Cholesky factor'
                         .Call(merPredDRXdiag, ptr())
                     },
                     b            = function(fac) {
                         'random effects on original scale for step factor fac'
                         .Call(merPredDb, ptr(), as.numeric(fac))
                     },
                     beta         = function(fac) {
                         'fixed-effects coefficients for step factor fac'
                         .Call(merPredDbeta, ptr(), as.numeric(fac))
                     },
                     copy         = function(shallow = FALSE) {
                         def <- .refClassDef
                         selfEnv <- as.environment(.self)
                         vEnv    <- new.env(parent=emptyenv())
                         
                         for (field in setdiff(names(def@fieldClasses), "Ptr")) {
                             if (shallow)
                                 assign(field, get(field, envir = selfEnv), envir = vEnv)
                             else {
                                 current <- get(field, envir = selfEnv)
                                 if (is(current, "envRefClass"))
                                     current <- current$copy(FALSE)
                                 ## hack (https://stat.ethz.ch/pipermail/r-devel/2014-March/068448.html)
                                 ## ... to ensure real copying
                                 ## forceCopy() does **NOT** work here, but +0 does
                                 ## we can get away with this because all fields other than Ptr are numeric
                                 assign(field, current+0, envir = vEnv)
                             }
                         }
                         do.call(merPredD$new, c(as.list(vEnv), n=nrow(vEnv$V), Class=def))
                     },
                     ldL2         = function() {
                         'twice the log determinant of the sparse Cholesky factor'
                         .Call(merPredDldL2, ptr())
                     },
                     ldRX2        = function() {
                         'twice the log determinant of the downdated dense Cholesky factor'
                         .Call(merPredDldRX2, ptr())
                     },
                     unsc         = function() {
                         'the unscaled variance-covariance matrix of the fixed-effects parameters'
                         .Call(merPredDunsc, ptr())
                     },
                     linPred      = function(fac) {
                         'evaluate the linear predictor for step factor fac'
                         .Call(merPredDlinPred, ptr(), as.numeric(fac))
                     },
                     installPars  = function(fac) {
                         'update u0 and beta0 to the values for step factor fac'
                         .Call(merPredDinstallPars, ptr(), as.numeric(fac))
                     },
                    initializePtr = function() {
                        Ptr <<- .Call(merPredDCreate, as(X, "matrix"), Lambdat,
                                      LamtUt, Lind, RZX, Ut, Utr, V, VtV, Vtr,
                                      Xwts, Zt, beta0, delb, delu, theta, u0)
                        .Call(merPredDsetTheta, Ptr, theta)
                        .Call(merPredDupdateXwts, Ptr, Xwts)
                        .Call(merPredDupdateDecomp, Ptr, NULL)
                     },
                     ptr          = function() {
                         'returns the external pointer, regenerating if necessary'
                         if (length(theta)) {
                             if (.Call(isNullExtPtr, Ptr)) initializePtr()
                         }
                         Ptr
                     },
                     setBeta0     = function(beta0) {
                         'install a new value of theta'
                         .Call(merPredDsetBeta0, ptr(), as.numeric(beta0))
                     },
                     setTheta     = function(theta) {
                         'install a new value of theta'
                         .Call(merPredDsetTheta, ptr(), as.numeric(theta))
                     },
                     solve        = function() {
                         'solve for the coefficient increments delu and delb'
                         .Call(merPredDsolve, ptr())
                     },
                     solveU       = function() {
                         'solve for the coefficient increment delu only (beta is fixed)'
                         .Call(merPredDsolveU, ptr())
                     },
                     setDelu      = function(val) {
                         'set the coefficient increment delu'
                         .Call(merPredDsetDelu , ptr(), as.numeric(val))
                     },
                     setDelb      = function(val) {
                         'set the coefficient increment delb'
                         .Call(merPredDsetDelb , ptr(), as.numeric(val))
                     },
                     sqrL         = function(fac) {
                         'squared length of u0 + fac * delu'
                         .Call(merPredDsqrL, ptr(), as.numeric(fac))
                     },
                     u            = function(fac) {
                         'orthogonal random effects for step factor fac'
                         .Call(merPredDu, ptr(), as.numeric(fac))
                     },
                     updateDecomp = function(XPenalty = NULL) {
                         'update L, RZX and RX from Ut, Vt and VtV'
                         invisible(.Call(merPredDupdateDecomp, ptr(), XPenalty))
                     },
                     updateL = function() {
                         'update LamtUt and L'
                         .Call(merPredDupdateL, ptr())
                     },
                     updateLamtUt = function() {
                         'update LamtUt and L'
                         .Call(merPredDupdateLamtUt, ptr())
                     },
                     updateRes    = function(wtres) {
                         'update Vtr and Utr using the vector of weighted residuals'
                         .Call(merPredDupdateRes, ptr(), as.numeric(wtres))
                     },
                     updateXwts   = function(wts) {
                         'update Ut and V from Zt and X using X weights'
                         .Call(merPredDupdateXwts, ptr(), wts)
                     }
                     )
                )
merPredD$lock("Lambdat", "LamtUt", "Lind", "RZX", "Ut", "Utr", "V", "VtV", "Vtr",
              "X", "Xwts", "Zt", "beta0", "delb", "delu", "theta", "u0")

##' Class \code{"merPredD"} - a dense predictor reference class
##'
##' A reference class for a mixed-effects model predictor module with a dense
##' model matrix for the fixed-effects parameters.  The reference class is
##' associated with a C++ class of the same name.  As is customary, the
##' generator object, \code{\link{merPredD}}, for the class has the same name as
##' the class.
##' @name merPredD-class
##' @note Objects from this reference class correspond to objects in a C++
##'     class.  Methods are invoked on the C++ class object using the external
##'     pointer in the \code{Ptr} field.  When saving such an object the external
##'     pointer is converted to a null pointer, which is why there are redundant
##'     fields containing enough information as R objects to be able to regenerate
##'     the C++ object.  The convention is that a field whose name begins with an
##'     upper-case letter is an R object and the corresponding field, whose name
##'     begins with the lower-case letter is a method.  References to the
##'     external pointer should be through the method, not directly through the
##'     \code{Ptr} field.
##' @section Extends: All reference classes extend and inherit methods from
##'     \code{"\linkS4class{envRefClass}"}.
##' @seealso \code{\link{lmer}}, \code{\link{glmer}}, \code{\link{nlmer}},
##'     \code{\link{merPredD}}, \code{\linkS4class{merMod}}.
##' @keywords classes
##' @examples
##'
##' showClass("merPredD")
##' str(slot(lmer(Yield ~ 1|Batch, Dyestuff), "pp"))
##'
NULL

##' Generator objects for the response classes
##'
##' The generator objects for the \code{\linkS4class{lmResp}},
##' \code{\linkS4class{lmerResp}}, \code{\linkS4class{glmResp}} and
##' \code{\linkS4class{nlsResp}} reference classes. Such objects are
##' primarily used through their \code{new} methods.
##'
##' @aliases lmResp lmerResp glmResp nlsResp
##' @param ... List of arguments (see Note).
##' @note Arguments to the \code{new} methods must be named arguments.
##' \itemize{
##' \item{y}{ the numeric response vector}
##' \item{family}{ a \code{\link{family}} object}
##' \item{nlmod}{ the nonlinear model function}
##' \item{nlenv}{ an environment holding data objects for evaluation of \code{nlmod}}
##' \item{pnames}{ a character vector of parameter names}
##' \item{gam}{ a numeric vector - the initial linear predictor}
##'}
##' @section Methods:
##' \describe{
##'     \item{\code{new(y=y)}:}{Create a new
##'          \code{\linkS4class{lmResp}} or \code{\linkS4class{lmerResp}} object.}
##'     \item{\code{new(family=family, y=y)}:}{Create a new
##'          \code{\linkS4class{glmResp}} object.}
##'     \item{\code{new(y=y, nlmod=nlmod, nlenv=nlenv, pnames=pnames,
##'                     gam=gam)}:}{Create a new
##'           \code{\linkS4class{nlsResp}} object.}
##' }
##' @seealso \code{\linkS4class{lmResp}}, \code{\linkS4class{lmerResp}},
##' \code{\linkS4class{glmResp}}, \code{\linkS4class{nlsResp}}
##' @keywords classes
##' @export
lmResp <-                               # base class for response modules
    setRefClass("lmResp",
                fields =
                list(Ptr     = "externalptr",
                     mu      = "numeric",
                     offset  = "numeric",
                     sqrtXwt = "numeric",
                     sqrtrwt = "numeric",
                     weights = "numeric",
                     wtres   = "numeric",
                     y       = "numeric"),
                methods =
                list(
                    allInfo = function() {
                        'return all the information available on the object'
                        data.frame(y=y, offset=offset, weights=weights, mu=mu,
                                   rwt=sqrtrwt, wres=wtres, Xwt=sqrtXwt)
                    },
                    initialize = function(...) {
                        if (!nargs()) return()
                        ll <- list(...)
                        if (is.null(ll$y)) stop("y must be specified")
                        y <<- as.numeric(ll$y)
                        n <- length(y)
                        mu <<- if (!is.null(ll$mu))
                            as.numeric(ll$mu) else numeric(n)
                        offset <<- if (!is.null(ll$offset))
                            as.numeric(ll$offset) else numeric(n)
                        weights <<- if (!is.null(ll$weights))
                            as.numeric(ll$weights) else rep.int(1,n)
                        sqrtXwt <<- if (!is.null(ll$sqrtXwt))
                            as.numeric(ll$sqrtXwt) else sqrt(weights)
                        sqrtrwt <<- if (!is.null(ll$sqrtrwt))
                            as.numeric(ll$sqrtrwt) else sqrt(weights)
                        wtres   <<- sqrtrwt * (y - mu)
                    },
                    copy         = function(shallow = FALSE) {
                        def <- .refClassDef
                        selfEnv <- as.environment(.self)
                        vEnv    <- new.env(parent=emptyenv())
                        for (field in setdiff(names(def@fieldClasses), "Ptr")) {
                            if (shallow)
                                assign(field, get(field, envir = selfEnv), envir = vEnv)
                            else {
                                current <- get(field, envir = selfEnv)
                                if (is(current, "envRefClass"))
                                    current <- current$copy(FALSE)
                                ## deep-copy hack +0
                                assign(field, forceCopy(current), envir = vEnv)
                            }
                        }
                        do.call("new", c(as.list(vEnv), Class=def))
                    },
                    initializePtr = function() {
                        Ptr <<- .Call(lm_Create, y, weights, offset, mu, sqrtXwt,
                                      sqrtrwt, wtres)
                        .Call(lm_updateMu, Ptr, mu)
                    },
                    ptr       = function() {
                        'returns the external pointer, regenerating if necessary'
                        if (length(y)) {
                            if (.Call(isNullExtPtr, Ptr)) initializePtr()
                        }
                        Ptr
                    },
                    setOffset  = function(oo) {
                        'change the offset in the model (used in profiling)'
                        .Call(lm_setOffset, ptr(), as.numeric(oo))
                    },
                    setResp    = function(rr) {
                        'change the response in the model, usually after a deep copy'
                        .Call(lm_setResp, ptr(), as.numeric(rr))
                    },
                    setWeights = function(ww) {
                        'change the prior weights in the model'
                        .Call(lm_setWeights, ptr(), as.numeric(ww))
                    },
                    updateMu  = function(gamma) {
                        'update mu, wtres and wrss from the linear predictor'
                        .Call(lm_updateMu, ptr(), as.numeric(gamma))
                    },
                    wrss      = function() {
                        'returns the weighted residual sum of squares'
                        .Call(lm_wrss, ptr())
                    })
                )

lmResp$lock("mu", "offset", "sqrtXwt", "sqrtrwt", "weights", "wtres")#, "y")

##' Classes \code{"lmResp"}, \code{"glmResp"}, \code{"nlsResp"} and
##' \code{"lmerResp"}
##'
##' Reference classes for response modules, including linear models,
##' \code{"lmResp"}, generalized linear models, \code{"glmResp"}, nonlinear
##' models, \code{"nlsResp"} and linear mixed-effects models, \code{"lmerResp"}.
##' Each reference class is associated with a C++ class of the same name.  As is
##' customary, the generator object for each class has the same name as the
##' class.
##' @name lmResp-class
##' @aliases lmResp-class glmResp-class lmerResp-class nlsResp-class
##' @note Objects from these reference classes correspond to objects in C++
##'     classes.  Methods are invoked on the C++ classes using the external pointer
##'     in the \code{ptr} field.  When saving such an object the external pointer is
##'     converted to a null pointer, which is why there are redundant fields
##'     containing enough information as R objects to be able to regenerate the C++
##'     object.  The convention is that a field whose name begins with an upper-case
##'     letter is an R object and the corresponding field whose name begins with the
##'     lower-case letter is a method.  Access to the external pointer should be
##'     through the method, not through the field.
##' @section Extends: All reference classes extend and inherit methods from
##'     \code{"\linkS4class{envRefClass}"}.  Furthermore, \code{"glmResp"},
##'     \code{"nlsResp"} and \code{"lmerResp"} all extend the \code{"lmResp"} class.
##' @seealso \code{\link{lmer}}, \code{\link{glmer}}, \code{\link{nlmer}},
##'     \code{\linkS4class{merMod}}.
##' @examples
##'
##' showClass("lmResp")
##' str(lmResp$new(y=1:4))
##' showClass("glmResp")
##' str(glmResp$new(family=poisson(), y=1:4))
##' showClass("nlsResp")
##' showClass("lmerResp")
##' str(lmerResp$new(y=1:4))
##' @keywords classes
NULL

##' @export
lmerResp <-
    setRefClass("lmerResp",
                fields=
                list(REML="integer"),
                contains="lmResp",
                methods=
                list(initialize = function(...) {
                         REML <<- as.integer(list(...)$REML)
                         if (length(REML) != 1L) REML <<- 0L
                         callSuper(...)
                     },
                     initializePtr = function() {
                         Ptr <<- .Call(lmer_Create, y, weights, offset, mu, sqrtXwt,
                                       sqrtrwt, wtres)
                         .Call(lm_updateMu, Ptr, mu - offset)
                         .Call(lmer_setREML, Ptr, REML)
                     },
                     ptr        = function() {
                         'returns the external pointer, regenerating if necessary'
                         if (length(y))
                             if (.Call(isNullExtPtr, Ptr)) initializePtr()
                         Ptr
                     },
                     objective  = function(ldL2, ldRX2, sqrL, sigma.sq = NULL) {
                         'returns the profiled deviance or REML criterion'
                         .Call(lmer_Laplace, ptr(), ldL2, ldRX2, sqrL, sigma.sq)
                     })
                )

setOldClass("family")

##' @export
glmResp <-
    setRefClass("glmResp",
                fields=
                list(eta="numeric", family="family", n="numeric"),
                contains="lmResp",
                methods=
                list(initialize = function(...) {
                         callSuper(...)
                         ll <- list(...)
                         if (is.null(ll$family)) stop("family must be specified")
                         family <<- ll$family
                         n <<- if (!is.null(ll$n)) as.numeric(ll$n) else rep.int(1,length(y))
                         eta <<- numeric(length(y))
                     },
                     aic          = function() {
                         .Call(glm_aic, ptr())
                     },
                     allInfo = function() {
                         'return all the information available on the object'
                         cbind(callSuper(),
                               data.frame(eta=eta, muEta=muEta(), var=variance(),
                                          WrkWt=sqrtWrkWt(), wrkRes=wrkResids(),
                                          wrkResp=wrkResp(), devRes=devResid()))
                     },
                     devResid  = function() {
                         'returns the vector of deviance residuals'
                         .Call(glm_devResid, ptr())
                     },
                     fam       = function() {
                         'returns the name of the glm family'
                         .Call(glm_family, ptr())
                     },
                     Laplace   = function(ldL2, ldRX2, sqrL) {
                         'returns the Laplace approximation to the profiled deviance'
                         .Call(glm_Laplace, ptr(), ldL2, ldRX2, sqrL)
                     },
                     link      = function() {
                         'returns the name of the glm link'
                         .Call(glm_link, ptr())
                     },
                     muEta     = function() {
                         'returns the diagonal of the Jacobian matrix, d mu/d eta'
                         .Call(glm_muEta, ptr())
                     },
                     ptr       = function() {
                         'returns the external pointer, regenerating if necessary'
                         if (length(y)) {
                             if (.Call(isNullExtPtr, Ptr)) {
                                 Ptr <<- .Call(glm_Create, family, y, weights, offset, mu, sqrtXwt,
                                               sqrtrwt, wtres, eta, n)
                                 .Call(glm_updateMu, Ptr, eta - offset)
                             }
                         }
                         Ptr
                     },
                     resDev = function() {
                         'returns the sum of the deviance residuals'
                         .Call(glm_resDev, ptr())
                     },
                     setTheta = function(theta) {
                         'sets a new value of theta, for negative binomial distribution only'
                         .Call(glm_setTheta, ptr(), as.numeric(theta))
                     },
                     sqrtWrkWt = function() {
                         'returns the square root of the working X weights'
                         .Call(glm_sqrtWrkWt, ptr())
                     },
                     theta = function() {
                         'query the value of theta, for negative binomial distribution only'
                         .Call(glm_theta, ptr())
                     },
                     updateMu = function(gamma) {
                         'update mu, residuals, weights, etc. from the linear predictor'
                         .Call(glm_updateMu, ptr(), as.numeric(gamma))
                     },
                     updateWts = function() {
                         'update the residual and X weights from the current value of eta'
                         .Call(glm_updateWts, ptr())
                     },
                     variance = function() {
                         'returns the vector of variances'
                         .Call(glm_variance, ptr())
                     },
                     wtWrkResp = function() {
                         'returns the vector of weighted working responses'
                         .Call(glm_wtWrkResp, ptr())
                     },
                     wrkResids = function() {
                         'returns the vector of working residuals'
                         .Call(glm_wrkResids, ptr())
                     },
                     wrkResp = function() {
                         'returns the vector of working responses'
                         .Call(glm_wrkResp, ptr())
                     }
                     )
                )


glmResp$lock("family", "n", "eta")

##' @export
nlsResp <-
    setRefClass("nlsResp",
                fields=
                list(gam="numeric",
                     nlmod="formula",
                     nlenv="environment",
                     pnames="character"
                     ),
                contains="lmResp",
                methods=
                list(initialize = function(...) {
                         callSuper(...)
                         ll <- list(...)
                         if (is.null(ll$nlmod)) stop("nlmod must be specified")
                         nlmod <<- ll$nlmod
                         if (is.null(ll$nlenv)) stop("nlenv must be specified")
                         nlenv <<- ll$nlenv
                         if (is.null(ll$pnames)) stop("pnames must be specified")
                         pnames <<- ll$pnames
                         if (is.null(ll$gam)) stop("gam must be specified")
                         stopifnot(length(ll$gam) == length(offset))
                         gam <<- ll$gam
                     },
                     Laplace =function(ldL2, ldRX2, sqrL) {
                         'returns the profiled deviance or REML criterion'
                         .Call(nls_Laplace, ptr(), ldL2, ldRX2, sqrL)
                     },
                     ptr     = function() {
                         'returns the external pointer, regenerating if necessary'
                         if (length(y)) {
                             if (.Call(isNullExtPtr, Ptr)) {
                                 Ptr <<- .Call(nls_Create, y, weights, offset, mu, sqrtXwt,
                                               sqrtrwt, wtres, gam, nlmod[[2]], nlenv, pnames)
                                 .Call(nls_updateMu, Ptr, gam)
                             }
                         }
                         Ptr
                     },
                     updateMu=function(gamma) {
                         'update mu, residuals, gradient, etc. given the linear predictor matrix'
                         .Call(nls_updateMu, ptr(), as.numeric(gamma))
                     })
                )

nlsResp$lock("nlmod", "nlenv", "pnames")

##' Generator object for the \code{\linkS4class{glmFamily}} class
##'
##' The generator object for the \code{\linkS4class{glmFamily}} reference class.
##' Such an object is primarily used through its \code{new} method.
##'
##'
##' @param ... Named argument (see Note below)
##' @note Arguments to the \code{new} method must be named arguments.
##' @section Methods: \describe{
##'     \item{\code{new(family=family)}}{Create a new
##'        \code{\linkS4class{glmFamily}} object}
##' }
##' @seealso \code{\linkS4class{glmFamily}}
##' @keywords classes
##' @export
glmFamily <-                            # used in tests of family definitions
    setRefClass("glmFamily",
                fields=list(Ptr="externalptr", family="family"),
                methods=
                list(
                     aic = function(y, n, mu, wt, dev) {
                         'returns the value from the aic member function, which is actually the deviance'
                         nn <- length(y <- as.numeric(y))
                         stopifnot(length(n <- as.numeric(n)) == nn,
                                   length(mu <- as.numeric(mu)) == nn,
                                   length(wt <- as.numeric(wt)) == nn,
                                   all(wt >= 0),
                                   length(dev <- as.numeric(dev)) == 1L)
                         .Call(glmFamily_aic, ptr(), y, n, mu, wt, dev)
                     },
                     devResid = function(y, mu, wt) {
                         'applies the devResid function to y, mu and wt'
                         mu <- as.numeric(mu)
                         wt <- as.numeric(wt)
                         y  <- as.numeric(y)
                         stopifnot(length(mu) == length(wt),
                                   length(mu) == length(y),
                                   all(wt >= 0))
                         .Call(glmFamily_devResid, ptr(), y, mu, wt)
                     },
                     link = function(mu) {
                         'applies the (forward) link function to mu'
                         .Call(glmFamily_link, ptr(), as.numeric(mu))
                     },
                     linkInv = function(eta) {
                         'applies the inverse link function to eta'
                         .Call(glmFamily_linkInv, ptr(), as.numeric(eta))
                     },
                     muEta = function(eta) {
                         'applies the muEta function to eta'
                         .Call(glmFamily_muEta, ptr(), as.numeric(eta))
                     },
                     ptr   = function() {
                         if (length(family))
                             if (.Call(isNullExtPtr, Ptr))
                                 Ptr <<- .Call(glmFamily_Create, family)
                         Ptr
                     },
                     setTheta = function(theta) {
                         'sets a new value of theta, for negative binomial distribution only'
                         .Call(glmFamily_setTheta, ptr(), as.numeric(theta))
                     },
                     theta = function() {
                         'query the value of theta, for negative binomial distribution only'
                         .Call(glmFamily_theta, ptr())
                     },
                     variance = function(mu) {
                         'applies the variance function to mu'
                         .Call(glmFamily_variance, ptr(), as.numeric(mu))
                     })
                )
##' Class \code{"glmFamily"} - a reference class for \code{\link{family}}
##'
##' This class is a wrapper class for \code{\link{family}} objects specifying a
##' distibution family and link function for a generalized linear model
##' (\code{\link{glm}}).  The reference class contains an external pointer to a
##' C++ object representing the class.  For common families and link functions
##' the functions in the family are implemented in compiled code so they can be
##' accessed from other compiled code and for a speed boost.
##'
##'
##' @name glmFamily-class
##' @docType class
##' @note Objects from this reference class correspond to objects in a C++
##' class.  Methods are invoked on the C++ class using the external pointer in
##' the \code{Ptr} field.  When saving such an object the external pointer is
##' converted to a null pointer, which is why there is a redundant field
##' \code{ptr} that is an active-binding function returning the external
##' pointer.  If the \code{Ptr} field is a null pointer, the external pointer is
##' regenerated for the stored \code{family} field.
##' @section Extends: All reference classes extend and inherit methods from
##' \code{"\linkS4class{envRefClass}"}.
##' @seealso \code{\link{family}}, \code{\link{glmFamily}}
##' @keywords classes
##' @examples
##'
##' str(glmFamily$new(family=poisson()))
NULL

##' Generator object for the golden search optimizer class.
##'
##' The generator objects for the \code{\linkS4class{golden}} class of a scalar
##' optimizer for a parameter within an interval.  The optimizer uses reverse
##' communications.
##'
##' @param \dots additional, optional arguments.  None are used at present.
##' @note Arguments to the \code{new} methods must be named arguments.
##' \code{lower} and \code{upper} are the bounds for the scalar parameter; they must be finite.
##' @section Methods:
##' \describe{
##'      \item{\code{new(lower=lower, upper=upper)}}{Create a new
##'         \code{\linkS4class{golden}} object.}
##' }
##' @seealso \code{\linkS4class{golden}}
##' @keywords classes
##' @export
golden <-
    setRefClass("golden", # Reverse communication implementation of Golden Search
                fields =
                list(
                     Ptr     = "externalptr",
                     lowerbd = "numeric",
                     upperbd = "numeric"
                     ),
                methods =
                list(
                     initialize = function(lower, upper, ...) {
                         stopifnot(length(lower <- as.numeric(lower)) == 1L,
                                   length(upper <- as.numeric(upper)) == 1L,
                                   lower > -Inf,
                                   upper < Inf,
                                   lower < upper)
                         lowerbd <<- lower
                         upperbd <<- upper
                         Ptr <<- .Call(golden_Create, lower, upper)
                     },
                     ptr        =  function() {
                         if (length(lowerbd))
                             if (.Call(isNullExtPtr, Ptr))
                                 Ptr <<- .Call(golden_Create, lowerbd, upperbd)
                         Ptr
                     },
                     newf       = function(value) {
                         stopifnot(length(value <- as.numeric(value)) == 1L)
                         .Call(golden_newf, ptr(), value)
                     },
                     value      = function() .Call(golden_value, ptr()),
                     xeval      = function() .Call(golden_xeval, ptr()),
                     xpos       = function() .Call(golden_xpos, ptr())
                     )
            )
##' Class \code{"golden"}
##'
##' A reference class for a golden search scalar optimizer using reverse
##' communication.
##'
##'
##' @name golden-class
##' @docType class
##' @section Extends: All reference classes extend and inherit methods from
##'    \code{"\linkS4class{envRefClass}"}.
##' @keywords classes
##' @examples
##'
##' showClass("golden")
##'
NULL

## Generator object for the Nelder-Mead optimizer class  "NelderMead"
##
## A reference class for a Nelder-Mead simplex optimizer allowing box
## constraints on the parameters and using reverse communication.
NelderMead <-
    setRefClass("NelderMead", # Reverse communication implementation of Nelder-Mead simplex optimizer
                fields =
                list(
                     Ptr     = "externalptr",
                     lowerbd = "numeric",
                     upperbd = "numeric",
                     xstep   = "numeric",
                     xtol    = "numeric"
                     ),
                methods =
                list(
                     initialize = function(lower, upper, xst, x0, xt, ...) {
                         stopifnot((n <- length(lower <- as.numeric(lower))) > 0L,
                                   length(upper <- as.numeric(upper)) == n,
                                   all(lower < upper),
                                   length(xst <- as.numeric(xst)) == n,
                                   all(xst != 0),
                                   length(x0 <- as.numeric(x0)) == n,
                                   all(x0 >= lower),
                                   all(x0 <= upper),
                                   all(is.finite(x0)),
                                   length(xt <- as.numeric(xt)) == n,
                                   all(xt > 0))
                         lowerbd <<- lower
                         upperbd <<- upper
                         xstep   <<- xst
                         xtol    <<- xt
                         Ptr <<- .Call(NelderMead_Create, lowerbd, upperbd, xstep, x0, xtol)
                     },
                     ptr          = function() {
                         if (length(lowerbd))
                             if (.Call(isNullExtPtr, Ptr))
                                 Ptr <<- .Call(NelderMead_Create, lowerbd, upperbd, xstep, x0, xtol)
                         Ptr
                     },
                     newf         = function(value) {
                         stopifnot(length(value <- as.numeric(value)) == 1L)
                         .Call(NelderMead_newf, ptr(), value)
                     },
                     setForceStop = function(stp=TRUE) .Call(NelderMead_setForce_stop, ptr(), stp),
                     setFtolAbs   = function(fta)      .Call(NelderMead_setFtol_abs, ptr(), fta),
                     setFtolRel   = function(ftr)      .Call(NelderMead_setFtol_rel, ptr(), ftr),
                     setMaxeval   = function(mxev)     .Call(NelderMead_setMaxeval, ptr(), mxev),
                     setMinfMax   = function(minf)     .Call(NelderMead_setMinf_max, ptr(), minf),
                     setIprint    = function(iprint)   .Call(NelderMead_setIprint, ptr(), iprint),
                     value        = function()         .Call(NelderMead_value, ptr()),
                     xeval        = function()         .Call(NelderMead_xeval, ptr()),
                     xpos         = function()         .Call(NelderMead_xpos, ptr())
                     )
            )


##' Class "merMod" of Fitted Mixed-Effect Models
##'
##' A mixed-effects model is represented as a \code{\linkS4class{merPredD}} object
##' and a response module of a class that inherits from class
##' \code{\linkS4class{lmResp}}.  A model with a \code{\linkS4class{lmerResp}}
##' response has class \code{lmerMod}; a \code{\linkS4class{glmResp}} response
##' has class \code{glmerMod}; and a \code{\linkS4class{nlsResp}} response has
##' class \code{nlmerMod}.
##'
##' @name merMod-class
##' @aliases merMod-class lmerMod-class glmerMod-class nlmerMod-class merMod
##' show,merMod-method
##' anova.merMod coef.merMod deviance.merMod
##' fitted.merMod formula.merMod logLik.merMod
##' model.frame.merMod model.matrix.merMod print.merMod
##' show.merMod summary.merMod
##' terms.merMod update.merMod
##' vcov.merMod print.summary.merMod show.summary.merMod
##' summary.summary.merMod vcov.summary.merMod
##' @docType class
##' @section Objects from the Class: Objects are created by calls to
##' \code{\link{lmer}}, \code{\link{glmer}} or \code{\link{nlmer}}.
##' @seealso \code{\link{lmer}}, \code{\link{glmer}}, \code{\link{nlmer}},
##' \code{\linkS4class{merPredD}}, \code{\linkS4class{lmerResp}},
##' \code{\linkS4class{glmResp}}, \code{\linkS4class{nlsResp}}
##' @keywords classes
##' @examples
##'
##' showClass("merMod")
##' methods(class="merMod")
##' @export
setClass("merMod",
         representation(Gp      = "integer",
                        call    = "call",
                        frame   = "data.frame", # "model.frame" is not S4-ized yet
                        flist   = "list",
                        cnms    = "list",
                        lower   = "numeric",
                        theta   = "numeric",
                        beta    = "numeric",
                        u       = "numeric",
                        devcomp = "list",
                        pp      = "merPredD",
                        optinfo = "list"))

##' @export
setClass("lmerMod", representation(resp="lmerResp"), contains="merMod")

##' @export
setClass("glmerMod", representation(resp="glmResp"), contains="merMod")

##' @export
setClass("nlmerMod", representation(resp="nlsResp"), contains="merMod")

##' Generator object for the rePos (random-effects positions) class
##'
##' The generator object for the \code{\linkS4class{rePos}} class used
##' to determine the positions and orders of random effects associated
##' with particular random-effects terms in the model.
##' @param \dots Argument list (see Note).
##' @note Arguments to the \code{new} methods must be named arguments.
##' \code{mer}, an object of class \code{"\linkS4class{merMod}"}, is
##' the only required/expected argument.
##' @section Methods:
##' \describe{
##'      \item{\code{new(mer=mer)}}{Create a new
##'         \code{\linkS4class{rePos}} object.}
##' }
##' @seealso \code{\linkS4class{rePos}}
##' @keywords classes
##' @export
rePos <-
    setRefClass("rePos",
                fields =
                list(
                     cnms    = "list", # component names (components are terms within a RE term)
                     flist   = "list", # list of grouping factors used in the random-effects terms
                     ncols   = "integer", # number of components for each RE term
                     nctot   = "integer", # total number of components per factor
                     nlevs   = "integer", # number of levels for each unique factor
                     offsets = "integer", # points to where each term starts
                     terms   = "list" # list with one element per factor, indicating corresponding term
                     ),
                methods =
                list(
                     initialize = function(mer, ...) {
                         ##' asgn indicates unique elements of flist
                         ##' 
                         stopifnot((ntrms <- length(Cnms <- mer@cnms)) > 0L,
                                   (length(Flist <- mer@flist)) > 0L,
                                   length(asgn  <- as.integer(attr(Flist, "assign"))) == ntrms)
                         cnms    <<- Cnms
                         flist   <<- Flist
                         ncols   <<- unname(vapply(cnms, length, 0L))
                         nctot   <<- unname(as.vector(tapply(ncols, asgn, sum)))
                         nlevs   <<- unname(vapply(flist, function(el) length(levels(el)), 0L))
                         # why not replace the sapply with ncols*nlevs[asgn] ??
                         offsets <<- c(0L, cumsum(sapply(seq_along(asgn),
                                                         function(i) ncols[i] * nlevs[asgn[i]])))
                         terms   <<- lapply(seq_along(flist), function(i) which(asgn == i))
                     }
                     )
            )
##' Class \code{"rePos"}
##'
##' A reference class for determining the positions in the random-effects vector
##' that correspond to particular random-effects terms in the model formula
##'
##' @name rePos-class
##' @docType class
##' @section Extends: All reference classes extend and inherit methods from
##'    \code{"\linkS4class{envRefClass}"}.
##' @keywords classes
##' @examples
##'
##' showClass("rePos")
##'
rePos$lock("cnms", "flist", "ncols", "nctot", "nlevs", "terms")

vcRep <-
    setRefClass("vcRep",
                fields =
                list(
                     theta     = "numeric",
                     lower     = "numeric",
                     Lambdat   = "dgCMatrix",
                     Lind      = "integer",
                     Gp        = "integer",
                     flist     = "list",
                     cnms      = "list",
                     ncols     = "integer",
                     nctot     = "integer",
                     nlevs     = "integer",
                     offsets   = "integer",
                     terms     = "list",
                     sig       = "numeric",
                     nms       = "character",
                     covar     = "list",
                     useSc     = "logical"
                     ),
                methods =
                list(
                     initialize = function(mer, ...) {
                         stopifnot((ntrms <- length(Cnms <- mer@cnms)) > 0L,
                                   (length(Flist <- mer@flist)) > 0L,
                                   length(asgn  <- as.integer(attr(Flist, "assign"))) == ntrms)
                         lower   <<- getME(mer, "lower")
                         theta   <<- getME(mer, "theta")
                         Lambdat <<- getME(mer, "Lambdat")
                         Lind    <<- getME(mer, "Lind")
                         Gp      <<- getME(mer, "Gp")
                         cnms    <<- Cnms
                         flist   <<- Flist
                         ncols   <<- unname(vapply(cnms, length, 0L))
                         nctot   <<- unname(as.vector(tapply(ncols, asgn, sum)))
                         nlevs   <<- unname(vapply(flist, function(el) length(levels(el)), 0L))
                         offsets <<- c(0L, cumsum(sapply(seq_along(asgn),
                                                         function(i) ncols[i] * nlevs[asgn[i]])))
                         terms   <<- lapply(seq_along(Flist), function(i) which(asgn == i))
                         sig     <<- sigma(mer)
                         nms     <<- names(Flist)[asgn]
                         covar   <<- mkVarCorr(sig, cnms, ncols, theta, nms)
                         useSc   <<- as.logical(getME(mer, "devcomp")$dims['useSc'])
                     },
                     asCovar     = function() {
                         ans <- lapply(covar,
                                       function(x) {
                                           attr(x, "correlation") <- attr(x, "stddev") <- NULL
                                           x
                                       })
                         attr(ans, "residVar") <- attr(covar, "sc")^2
                         ans
                     },
                     asCorr      = function() {
                         ans <- lapply(covar, function(x)
                                       list(correlation=attr(x, "correlation"),
                                            stddev=attr(x, "stddev")))
                         attr(ans, "residSD") <- attr(covar, "sc")
                         ans
                     },
                     setTheta    = function(ntheta) {
                         stopifnot(length(ntheta <- as.numeric(ntheta)) == length(lower),
                                   all(ntheta >= lower))
                         theta   <<- ntheta
                         covar   <<- mkVarCorr(sig, cnms, ncols, theta, nms)
                     },
                     setSc       = function(nSc) {
                         stopifnot(useSc,
                                   length(nSc <- as.numeric(nSc)) == 1L)
                         sig     <<- nSc
                         covar   <<- mkVarCorr(sig, cnms, ncols, theta, nms)
                     },
                     setResidVar = function(nVar) setSc(sqrt(as.numeric(nVar))),
                     setRECovar  = function(CV) {
                         if (is.matrix(CV) && length(covar) == 1L) {
                             CV <- list(CV)
                             names(CV) <- names(covar)
                         }
                         covsiz <- sapply(covar, ncol)
                         stopifnot(is.list(CV),
                                   all(names(CV) == names(covar)),
                                   all(sapply(CV, isSymmetric)),
                                   all(sapply(CV, ncol) == covsiz))
                         if (!all(vapply(cnms, length, 1L) == covsiz))
                             error("setRECovar currently requires distinct grouping factors")
                         theta <<- sapply(CV, function(mm)
                                      {
                                          ff <- t(chol(mm))/sig
                                          ff[upper.tri(ff, diag=TRUE)]
                                      })
                     })
                )
