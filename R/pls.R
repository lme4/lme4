##' Compute linear mixed model deviance
##'
##' A pure \code{R} implementation of Doug Bates'
##' penalized least squares (PLS) approach for computing
##' linear mixed model deviances. The purpose
##' is to clarify how PLS works without having
##' to read through C++ code, and as a sandbox for
##' trying out modifications to PLS.
##'
##' @param obj output of \code{lFormula} or a model matrix
##' @param theta covariance parameters
##' @param y response
##' @param weights prior weights
##' @param offset offset
##' @param REML calculate REML deviance?
##'
##' @return deviance with a bunch of attributes:
##' \itemize{
##' \item[devML] ML deviance
##' \item[devREML] REML deviance
##' \item[beta] fixed effects at convergence
##' \item[u] conditional mode at convergence
##' \item[sigmaML] ML estimate of the residual standard deviation 
##' \item[sigmaREML] REML estimate of the residual standard deviation
##' }
##' @export
pls <- function(obj,theta,y,...) UseMethod("pls")

##' @rdname pls
##' @method pls matrix
##' @S3method pls matrix
pls.matrix <- function(obj,theta,y,Zt,Lambdat,Lind,
                       weights = rep(1, length(y)),
                       offset = rep(0, length(y)),
                       REML = TRUE){
    n <- nrow(obj)
    p <- ncol(obj)
    q <- nrow(Zt)
                                        # update Lambdat
    Lambdat@x <- theta[Lind]
                                        # calculate weighted design matrices
    W <- Diagonal(x = sqrt(weights))
    Ut <- Zt %*% W
    LamtUt <- Lambdat %*% Ut
    V <- W %*% obj
                                        # calculate weighted response
    wtResp <- as.vector(W %*% y)
    ## sparse Cholesky decomposition of cross-product matrix
    L <- Cholesky(tcrossprod(LamtUt), perm=FALSE, LDL = FALSE, Imult=1)
    RZX <- solve(L, LamtUt%*%V, system = "L")
    RX <- try(chol(crossprod(V) - crossprod(RZX)))
    if(inherits(RX, "try-error")) return(NA)
                                        # solve for fixed effect estimates and conditional modes
    cu <- solve(L, LamtUt %*% wtResp, system = "L")
    cb <- solve(t(RX), crossprod(V,wtResp) - crossprod(RZX,cu))
    beta <- solve(RX, cb)
    u <- solve(L, cu - RZX %*% beta, system = "Lt")
                                        # calculate conditional mean of the response
    mu <- crossprod(Zt,crossprod(Lambdat,u)) + (obj %*% beta) + offset
                                        # calculate weighted residuals
    wtres <- sqrt(weights)*(y-mu)
                                        # calculate residual sums of squares
    wrss <- sum(wtres^2)
    pwrss <- wrss + sum(u^2)            # penalize
                                        # calculate log determinants
    ldL2 <- 2*determinant(L,logarithm=TRUE)$modulus
    ldRX2 <- 2*determinant(RX,logarithm=TRUE)$modulus
    attributes(ldL2) <- attributes(ldRX2) <- NULL
                                        # calculate profiled deviance for theta
    devML <- ldL2 + n*(1 + log(2*pi*pwrss) - log(n))
    devREML <- ldL2 + ldRX2 + (n-p)*(1 + log(2*pi*pwrss) -
                                     log(n-p))
    if(REML) dev <- devREML
    else dev <- devML
    sigmaML <- sqrt(pwrss/sum(weights))
    return(structure(dev, devML = devML, devREML = devREML, beta = beta, u = u,
                     sigmaML = sigmaML,
                     sigmaREML = sigmaML*(n/(n-p))))
}
##' @rdname pls
##' @method pls list
##' @S3method pls list
pls.list <- function(obj,theta,y,weights=rep(1,length(y)),
                     offset=rep(0,length(y)),REML=TRUE) {
    retrm <- obj$reTrms
    pls.matrix(obj$X,theta,y,retrm$Zt,retrm$Lambdat,retrm$Lind,weights,offset,REML)
}
