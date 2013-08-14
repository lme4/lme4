##' Penalized iteratively reweighted least squares
##'
##' A pure \code{R} implementation of Doug Bates'
##' penalized iteratively reweighted least squares (PIRLS)
##' algorithm for computing generalized linear mixed model
##' deviances. The purpose is to clarify how
##' PIRLS works without having to read through C++ code,
##' and as a sandbox for trying out modified versions of
##' PIRLS.
##'
##' @param theta covariance parameters
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
pirls <- function(theta, beta, u,
                  mu, eta,
                  glmod, y,
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

 # initialize beta and u
 # (only necessary if step-halving starts immediately)
 if(missing(beta)) beta <- rep(0, p)
 if(missing(u)) u <- rep(0, q)
 newbeta <- beta
 
 # update theta
 Lambdat@x <- theta[Lind]

 # initialize penalized deviance and convergence flag
 oldpdev <- .Machine$double.xmax
 cvgd <- FALSE

 # PIRLS
 for(i in 1:npirls){
 
  # update w and muEta
  w <- weights/as.vector(family()$variance(mu))
  muEta <- as.vector(family()$mu.eta(eta))

  # update Ut and V
  Xwts <- sqrt(w)*muEta
  Ut <- Lambdat%*%Zt%*%Diagonal(n, Xwts)
  UtU <- tcrossprod(Ut)
  if(!nAGQ){
    V <- diag(Xwts, n, n)%*%X
    VtV <- t(V)%*%V
    UtV <- Ut%*%V
  }

  # update L, RZX, and RX
  L <- Cholesky(UtU, LDL = FALSE, Imult=1)
  if(!nAGQ){
    RZX <- solve(L, UtV, system = "L")
    RX <- try(chol(VtV - crossprod(RZX)))
    if(inherits(RX, "try-error"))
      stop("Downdated VtV not positive definite")
  }

  # update Utr and Vtr
  if(nAGQ){
    r <- Xwts*(eta - offset - X%*%beta + ((y-mu)/muEta))    
  }
  else {
    r <- Xwts*(eta - offset + ((y-mu)/muEta))
    Vtr <- t(V)%*%r
  }
  Utr <- Ut%*%r

  # solve for u and beta
  cu <- solve(L, Utr, system = "L")
  if(!nAGQ){
    cb <- solve(t(RX), Vtr - t(RZX)%*%cu)
    newbeta <- as.vector(solve(RX, cb))
    newu <- as.vector(solve(L, cu - RZX%*%newbeta,
                            system = "Lt"))
  }
  else {
    newu <- as.vector(solve(L, cu, system = "Lt"))
  }

  # update mu and eta
  eta <- offset + as.vector((t(Lambdat%*%Zt)%*%newu) +
                            (X%*%newbeta))
  mu <- as.vector(family()$linkinv(eta))
  
  # compute penalized deviance
  pdev <- sum(family()$dev.resid(y, mu, weights)) +
    sum(newu^2)
  if(abs((oldpdev - pdev) / pdev) < tol){
    cvgd <- TRUE
    break
  }

  # step-halving
  if(pdev > oldpdev){
    for(j in 1:10){
      newu <- (newu + u)/2
      if(!nAGQ) newbeta <- (newbeta + beta)/2
      eta <- offset +
        as.vector((t(Lambdat%*%Zt)%*%newu) +
                  (X%*%newbeta))
      mu <- as.vector(family()$linkinv(eta))
      pdev <- sum(family()$dev.resid(y, mu, weights)) +
        sum(newu^2)
      if(!(pdev > oldpdev)) break
    }
    if((pdev - oldpdev) > tol)
      stop("Step-halving failed")
  }
  oldpdev <- pdev
  u <- newu
  if(!nAGQ) beta <- newbeta
 }
 if(cvgd){
    # calculate laplace deviance approximation
    ldL2 <- 2*determinant(L, logarithm = TRUE)$modulus
    attributes(ldL2) <- NULL
    ldev <- pdev + ldL2 + (q/2)*log(2*pi)
    out <- structure(ldev, pdev = pdev,
                     beta = newbeta, u = newu)
    return(out)
 }
 stop("PIRLS failed to converge")
}

##' @rdname pirls
##' @export
pirls1 <- function(thetabeta, ...){
  p <- ncol(glmod$X)
  nth <- length(thetabeta) - p
  theta <- thetabeta[1:nth]
  beta <- thetabeta[-(1:nth)]
  pirls(theta, beta, ..., nAGQ = 1)
}
