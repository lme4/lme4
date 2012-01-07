### Class definitions for the package

setClass("lmList",
         representation(call = "call",
                        pool = "logical"),
         contains = "list")

## TODO: export
setClass("lmList.confint", contains = "array")

### FIXME
### shouldn't we have "merPred"  with two *sub* classes "merPredD" and "merPredS"
### for the dense and sparse X cases ?

merPredD <- 
    setRefClass("merPredD", # Predictor class for mixed-effects models with dense X
                fields =
                list(Lambdat = "dgCMatrix",
                     LamtUt  = "dgCMatrix",
                     Lind    = "integer",
                     Ptr     = "externalptr",
                     RZX     = "matrix",
                     Ut      = "dgCMatrix",
                     Utr     = "numeric",
                     V       = "matrix",
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
### FIXME: probably don't need S as Xwts is required for nlmer 
                     initialize = function(X, Zt, Lambdat, Lind, theta, S, ...) {
                         if (!nargs()) return
                         X <<- as(X, "matrix")
                         Zt <<- as(Zt, "dgCMatrix")
                         Lambdat <<- as(Lambdat, "dgCMatrix")
                         Lind <<- as.integer(Lind)
                         theta <<- as.numeric(theta)
                         N <- nrow(X)
                         S <- as.integer(S)[1]
                         n <- N %/% S
                         p <- ncol(X)
                         q <- nrow(Zt)
                         stopifnot(length(theta) > 0L,
                                   length(Lind) > 0L, 
                                   all(sort(unique(Lind)) == seq_along(theta)),
                                   S > 0L,
                                   n * S == N)
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
                         Ut <<- if (S == 1) Zt + 0 else
                             Zt %*% sparseMatrix(i=seq_len(N), j=as.integer(gl(n, 1, N)), x=rep.int(1,N))
                         LamtUt <<- Lambdat %*% Ut   
                         Xw <- list(...)$Xwts
                         Xwts <<- if (is.null(Xw)) rep.int(1, N) else as.numeric(Xw)
                         updateXwts(Xwts)
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
                     ptr          = function() {
                         'returns the external pointer, regenerating if necessary'
                         if (length(theta)) {
                             if (.Call(isNullExtPtr, Ptr)) {
                                 Ptr <<- .Call(merPredDCreate, as(X, "matrix"), Lambdat,
                                               LamtUt, Lind, RZX, Ut, Utr, V, VtV, Vtr,
                                               Xwts, Zt, beta0, delb, delu, theta, u0)
                                 .Call(merPredDsetTheta, Ptr, theta)
                                 .Call(merPredDupdateXwts, Ptr, Xwts)
                                 .Call(merPredDupdateDecomp, Ptr)
                             }
                         }
                         Ptr
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
                     sqrL         = function(fac) {
                         'squared length of u0 + fac * delu'
                         .Call(merPredDsqrL, ptr(), as.numeric(fac))
                     },
                     u            = function(fac) {
                         'orthogonal random effects for step factor fac'
                         .Call(merPredDu, ptr(), as.numeric(fac))
                     },
                     updateDecomp = function() {
                         'update L, RZX and RX from Ut, Vt and VtV'
                         .Call(merPredDupdateDecomp, ptr())
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
                     ptr       = function() {
                         'returns the external pointer, regenerating if necessary'
                         if (length(y)) {
                             if (.Call(isNullExtPtr, Ptr)) {
                                 Ptr <<- .Call(lm_Create, y, weights, offset, mu, sqrtXwt,
                                               sqrtrwt, wtres)
                                 .Call(lm_updateMu, Ptr, mu)
                             }
                         }
                         Ptr
                     },
                     setOffset = function(oo) {
                         'change the offset in the model (used in profiling)'
                         .Call(lm_setOffset, ptr(), as.numeric(oo))
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
                
lmResp$lock("mu", "offset", "sqrtXwt", "sqrtrwt", "weights", "wtres", "y")

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
                     ptr        = function() {
                         'returns the external pointer, regenerating if necessary'
                         if (length(y)) {
                             if (.Call(isNullExtPtr, Ptr)) {
                                 Ptr <<- .Call(lmer_Create, y, weights, offset, mu, sqrtXwt,
                                               sqrtrwt, wtres)
                                 .Call(lm_updateMu, Ptr, mu - offset)
                                 .Call(lmer_setREML, Ptr, REML)
                             }
                         }
                         Ptr
                     },
                     objective  = function(ldL2, ldRX2, sqrL) {
                         'returns the profiled deviance or REML criterion'
                         .Call(lmer_Laplace, ptr, ldL2, ldRX2, sqrL)
                     })
                )

setOldClass("family")

glmResp <-
    setRefClass("glmResp",
                fields=
                list(eta="numeric", family="family", n="integer"),
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
                     sqrtWrkWt = function() {
                         'returns the square root of the working X weights'
                         .Call(glm_sqrtWrkWt, ptr())
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
                     wrkResids = function() {
                         'returns the vector of working residuals'
                         .Call(glm_wrkResids, ptr())
                     },
                     wrkResp = function() {
                         'returns the vector of working residuals'
                         .Call(glm_wrkResp, ptr())
                     }
                     )
                )


glmResp$lock("family", "n", "eta")

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

glmFamily <-                            # used in tests of family definitions
    setRefClass("glmFamily",
                fields=list(Ptr="externalptr", family="family"),
                methods=
                list(
                     devResid = function(mu, weights, y) {
                         'applies the devResid function to mu, weights and y'
                         mu <- as.numeric(mu)
                         weights <- as.numeric(weights)
                         y <- as.numeric(y)
                         stopifnot(length(mu) == length(weights),
                                   length(mu) == length(y),
                                   all(weights >= 0))
                         .Call(glmFamily_devResid, ptr(), mu, weights, y)
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
                     variance = function(mu) {
                         'applies the variance function to mu'
                         .Call(glmFamily_variance, ptr(), as.numeric(mu))
                     })
                )

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
                     value        = function()         .Call(NelderMead_value, ptr()),
                     xeval        = function()         .Call(NelderMead_xeval, ptr()),
                     xpos         = function()         .Call(NelderMead_xpos, ptr())
                     )
            )

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
                        pp      = "merPredD"))

setClass("lmerMod", representation(resp="lmerResp"), contains="merMod")

setClass("glmerMod", representation(resp="glmResp"), contains="merMod")

setClass("nlmerMod", representation(resp="nlsResp"), contains="merMod")
