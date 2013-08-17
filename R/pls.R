##' @importMethodsFrom Matrix t %*% crossprod diag tcrossprod solve determinant
NULL

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
##' @param ... Arguments to pass to other functions
##'
##' @return a function that evaluates the deviance or REML criterion
##' @export
pls <- function(obj,y,...) UseMethod("pls")

##' @rdname pls
##' @param Zt transpose of the sparse model matrix for the random effects
##' @param Lambdat upper triangular sparse Cholesky factor of the
##'    relative covariance matrix of the random effects
##' @param thfun a function that takes a value of \code{theta} and produces
##'    the non-zero elements of \code{Lambdat}.  The structure of \code{Lambdat}
##'    cannot change, only the numerical values
##' @param weights prior weights
##' @param offset offset
##' @param REML calculate REML deviance?
##' @keywords models
##' @examples
##' library(lme4)
##' library(nloptwrap)
##' lmod <- lFormula(Reaction ~ Days + (Days|Subject), sleepstudy)
##' devf <- pls(lmod,sleepstudy$Reaction)
##' devf(c(1,0,1))             # starting value
##' bobyqa(c(1, 0, 1), devf, lower=c(0,-Inf,0))[c("par","value")]
##' @method pls matrix
##' @S3method pls matrix
pls.matrix <- function(obj,y,Zt,Lambdat,thfun,
                       weights = rep(1, n),
                       offset = numeric(n),
                       REML = TRUE,...){
    n <- length(y)
    p <- ncol(obj)
    q <- nrow(Zt)
    stopifnot(nrow(obj) == n,
              ncol(Zt) == n,
              nrow(Lambdat) == q,
              ncol(Lambdat) == q,
              is.function(thfun))
                                        # calculate weighted products
    W <- Diagonal(x = weights)
    L <- Cholesky(tcrossprod(Lambdat %*% Zt), LDL = FALSE, Imult=1)
    XtWX <- crossprod(obj, W %*% obj)
    XtWy <- crossprod(obj, W %*% y)
    ZtWX <- Zt %*% (W %*% obj)
    ZtWy <- Zt %*% (W %*% y)
    beta <- numeric(p)
    u <- numeric(q)
    mu <- numeric(n)
    DD <- NULL
    RZX <- matrix(0,nrow=q,ncol=p)
    cu <- numeric(q)
    function(theta) {
        Lambdat@x[] <<- thfun(theta) 
        L <- update(L, Lambdat %*% Zt, mult = 1)
        ## solve system from equation 30
        cu[] <- as.vector(solve(L, solve(L, Lambdat %*% ZtWy, system="P"),
                                system="L"))
        ## solve system from eqn. 31
        RZX[] <<- as.vector(solve(L, solve(L, Lambdat %*% ZtWX, system="P"),
                                  system="L"))
        ## downdate XtWX and form Cholesky factor (eqn. 32)
        DD <<- as(XtWX - crossprod(RZX), "dpoMatrix")
        ## conditional estimate of fixed-effects coefficients (solve eqn. 33)
        beta[] <<- as.vector(solve(DD, XtWy - crossprod(RZX, cu)))
        ## conditional mode of the spherical random-effects coefficients (eqn. 34)
        u[] <<- as.vector(solve(L, solve(L, cu - RZX %*% beta, system = "Lt"), system="Pt"))
                                        # conditional mean of the response

        ## crossprod(Zt,crossprod(Lambdat,u))  == Z Lambda u == Z b
        mu[] <<- as.vector(crossprod(Zt,crossprod(Lambdat,u)) + obj %*% beta + offset)
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
                     offset=rep(0,length(y)),REML=TRUE,...) {
    retrm <- obj$reTrms
    thfun <- retrm$thfun
    if (is.null(thfun)) thfun <- local({
        Lind <- retrm$Lind
        function(theta) theta[Lind]
    })
    pls.matrix(obj$X,y,retrm$Zt,retrm$Lambdat,thfun,weights,offset,REML)
}
