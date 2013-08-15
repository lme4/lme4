##' Create linear mixed model deviance function
##'
##' A pure \code{R} implementation of the
##' penalized least squares (PLS) approach for computing
##' linear mixed model deviances. The purpose
##' is to clarify how PLS works without having
##' to read through C++ code, and as a sandbox for
##' trying out modifications to PLS.
##'
##' @param obj output of \code{lFormula} or a model matrix
##' @param y response
##' @param weights prior weights
##' @param offset offset
##' @param REML calculate REML deviance?
##'
##' @return a function that evaluates the deviance or REML criterion
##' @export
pls <- function(obj,theta,y,...) UseMethod("pls")

##' @rdname pls
##' @method pls matrix
##' @S3method pls matrix
pls.matrix <- function(X,y,Zt,Lambdat,thfun,
                       weights = rep(1, length(y)),
                       offset = rep(0, length(y)),
                       REML = TRUE){
    n <- length(y)
    p <- ncol(X)
    q <- nrow(Zt)
    stopifnot(nrow(X) == n,
              ncol(Zt) == n,
              nrow(Lambdat) == q,
              ncol(Lambdat) == q,
              is.function(thfun))
                                        # calculate weighted products
    W <- Diagonal(x = weights)
    L <- Cholesky(tcrossprod(Lambdat %*% Zt), perm=FALSE, LDL = FALSE, Imult=1)
    XtWX <- crossprod(X, W %*% X)
    XtWy <- crossprod(X, W %*% y)
    ZtWy <- Zt %*% (W %*% y)
    ZtWX <- Zt %*% (W %*% X)
    beta <- numeric(ncol(X))
    u <- numeric(nrow(Zt))
    mu <- numeric(nrow(X))
    function(theta) {
        Lambdat@x[] <- thfun(theta) 
        L <- update(L, Lambdat %*% Zt, mult = 1)
        cu <- solve(L, Lambdat %*% ZtWy, system = "L")
        RZX <- solve(L, Lambdat %*% ZtWX, system = "L")
        ## downdate XtWX and form Cholesky factor
        DD <- as(XtWX - crossprod(RZX), "dpoMatrix")
        ## conditional estimate of fixed-effects coefficients
        beta <<- solve(DD, XtWy - crossprod(RZX, cu))
        ## conditional mode of the spherical random-effects coefficients
        u <<- solve(L, cu - RZX %*% beta, system = "Lt")
                                        # conditional mean of the response
        mu <<- crossprod(Zt,crossprod(Lambdat,u)) + X %*% beta + offset
                                        # weighted residuals
        wtres <- sqrt(weights)*(y-mu)
                                        # weighted residual sums of squares
        wrss <- sum(wtres^2)
        pwrss <- wrss + sum(u^2)        # penalize
                                        # log determinants
        ldL2 <- 2*determinant(L,logarithm=TRUE)$modulus
#        ldRX2 <- determinant(DD,logarithm=TRUE)$modulus # need method in Matrix
        attributes(ldL2) <- NULL #attributes(ldRX2) <- NULL
                                        # profiled deviance or REML criterion
        ## if (REML) 
        ##     ldL2 + ldRX2 + (n-p)*(1 + log(2*pi*pwrss) - log(n-p))
        ## else
        ldL2 + n*(1 + log(2*pi*pwrss) - log(n))
    }
}
##' @rdname pls
##' @method pls list
##' @S3method pls list
pls.list <- function(obj,y,weights=rep(1,length(y)),
                     offset=rep(0,length(y)),REML=TRUE) {
    retrm <- obj$reTrms
    thfun <- retrm$thfun
    if (is.null(thfun)) thfun <- local({
        Lind <- retrm$Lind
        function(theta) theta[Lind]
    })
    pls.matrix(obj$X,y,retrm$Zt,retrm$Lambdat,thfun,weights,offset,REML)
}
