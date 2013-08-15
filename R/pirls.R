##' Create an approximate deviance evaluation function for GLMMs
##'
##' A pure \code{R} implementation of the
##' penalized iteratively reweighted least squares (PIRLS)
##' algorithm for computing generalized linear mixed model
##' deviances. The purpose is to clarify how
##' PIRLS works without having to read through C++ code,
##' and as a sandbox for trying out modified versions of
##' PIRLS.
##'
##' @param beta initial beta
##' @param u initial u
##' @param mu fitted values
##' @param eta linear predictor
##' @param glmod output of \code{glFormula}
##' @param y response
##' @param family family
##' @param weights prior weights
##' @param offset offset
##' @param tol convergence tolerance
##' @param npirls maximum number of iterations
##' @param nAGQ either 0 (PIRLS for \code{u} and \code{beta}) or 1 (\code{u} only).
##'  currently no quadature is available
##' @param thetabeta \code{c(theta,beta)} (\code{pirls1} only)
##' @param ... Arguments to pass to \code{pirls}
##'
##' @details \code{pirls1} is a convenience function for optimizing \code{pirls}
##' under \code{nAGQ = 1}. In particular, it wraps \code{theta} and \code{beta}
##' into a single argument \code{thetabeta}.
##' 
##' @return the laplace approximated deviance with a bunch of attributes:
##' \itemize{
##' \item[pdev] penalized deviance at convergence
##' \item[beta] fixed effects at convergence (\code{pirls0} only)
##' \item[u] spherized random effects at convergence
##' }
##' @importMethodsFrom Matrix t %*% crossprod diag tcrossprod solve determinant
##' @rdname pirls
##' @export
pirls <- function(glmod, y, eta,
                  family = binomial,
                  weights = rep(1, length(y)),
                  offset = rep(0, length(y)),
                  tol = 10^-6, npirls = 30,
                  nAGQ = 0){
                                        # expand glmod
    Lind <- glmod$reTrms$Lind
    Lambdat <- glmod$reTrms$Lambdat
    Zt <- glmod$reTrms$Zt
    X <- glmod$X
    n <- nrow(X)
    p <- ncol(X)
    q <- nrow(Zt)
    if (is.function(family)) family <- family() # ensure family is a list
    mu <- family$linkinv(eta)
    beta <- numeric(p)
    u <- numeric(q)
    L <- Cholesky(tcrossprod(Lambdat %*% Zt), perm=FALSE, LDL=FALSE, Imult=1)
    betaind <- -seq_len(max(Lind))      # indices to drop 1:nth
    function(thetabeta) {
        Lambdat@x[] <<- thetabeta[Lind] # update Lambdat
        beta[] <<- thetabeta[betaind]
                                        # initialization
        oldpdev <- .Machine$double.xmax
        cvgd <- FALSE
        for(i in 1:npirls){             # PIRLS
                                        # update w and muEta
            W <- Diagonal(x=weights/family$variance(mu))
            muEta <- as.vector(family$mu.eta(eta))
                                        # update Ut and V
            Xwts <- Diagonal(x = sqrt(W@x) * muEta)
            Ut <- Lambdat %*% Zt %*% Xwts
            L <- update(L, Ut, 1)
            if(!nAGQ){                  # update Utr and Vtr
                V <- Xwts %*% Matrix(X)
                VtV <- crossprod(V)
                RZX <- solve(L, Ut %*% V, system = "L")
                DD <- as(as(VtV - crossprod(RZX), "denseMatrix"),"dpoMatrix")
                r <- Xwts %*% (eta - offset + ((y-mu)/muEta))
                cu <- solve(L, Ut %*% r, system = "L")
                newbeta <- as.vector(solve(DD, crossprod(V,r) - crossprod(RZX,cu)))
                newu <- as.vector(solve(L, cu - RZX %*% newbeta, system = "Lt"))
            } else {
                r <- Xwts %*% (eta - offset - X%*%beta + ((y-mu)/muEta))
                newu <- as.vector(solve(L, Ut %*% r))
                newbeta <- beta
            }
                                        # update mu and eta
            eta[] <<- as.vector(offset + X %*% newbeta +
                                crossprod(Zt,crossprod(Lambdat,newu)))
            mu[] <<- family$linkinv(eta)
                                        # compute penalized deviance
            pdev <- sum(family$dev.resid(y, mu, weights)) + sum(newu^2)
            if(abs((oldpdev - pdev) / pdev) < tol){
                cvgd <- TRUE
                break
            }
                                        # step-halving
            if(pdev > oldpdev){
                for(j in 1:10){
                    newu <- (newu + u)/2
                    if(!nAGQ) newbeta <- (newbeta + beta)/2
                    eta[] <- as.vector(offset + X %*% newbeta +
                                       crossprod(Zt, crossprod(Lambdat,newu)))
                    mu <- family$linkinv(eta)
                    pdev <- sum(family$dev.resid(y, mu, weights)) + sum(newu^2)
                    if(!(pdev > oldpdev)) break
                }
                if((pdev - oldpdev) > tol) stop("Step-halving failed")
            }
            oldpdev <- pdev
            u <- newu
            if(!nAGQ) beta[] <<- newbeta
        }
        if(!cvgd) stop("PIRLS failed to converge")
        ## calculate Laplace deviance approximation
        ldL2 <- 2*determinant(L, logarithm = TRUE)$modulus
        attributes(ldL2) <- NULL
        pdev + ldL2 + (q/2)*log(2*pi)
    }
}
