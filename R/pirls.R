##' Create an approximate deviance evaluation function for GLMMs using Laplace
##' Must use the flexLambda branch of lme4
##'
##' A pure \code{R} implementation of the
##' penalized iteratively reweighted least squares (PIRLS)
##' algorithm for computing generalized linear mixed model
##' deviances. The purpose is to clarify how
##' PIRLS works without having to read through C++ code,
##' and as a sandbox for trying out modified versions of
##' PIRLS.
##'
##' @param glmod output of \code{glFormula}
##' @param y response
##' @param eta linear predictor
##' @param family a \code{glm} family object
##' @param weights prior weights
##' @param offset offset
##' @param tol convergence tolerance
##' @param npirls maximum number of iterations
##' @param nAGQ either 0 (PIRLS for \code{u} and \code{beta}) or 1 (\code{u} only).
##'     currently no quadature is available
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
##' @export
pirls <- function(glmod, y, eta,
                  family = binomial,
                  weights = rep(1, length(y)),
                  offset = rep(0, length(y)),
                  tol = 10^-6, npirls = 30,
                  nAGQ = 1, verbose=0L){
                                        # expand glmod
    retrms <- glmod$reTrms
    thfun <- retrms$thfun
    nth <- max(thfun(retrms$theta))
    betaind <- -seq_len(nth) # indices to drop 1:nth
    Lambdat <- retrms$Lambdat
    Zt <- retrms$Zt
    X <- glmod$X
    n <- nrow(X)
    p <- ncol(X)
    q <- nrow(Zt)
    if (is.function(family)) family <- family() # ensure family is a list
    linkinv <- family$linkinv
    variance <- family$variance
    muEta <- family$mu.eta
    aic <- family$aic
    sqDevResid <- family$dev.resid
    mu <- linkinv(eta)
    beta <- numeric(p)
    u <- numeric(q)
    L <- Cholesky(tcrossprod(Lambdat %*% Zt), perm=FALSE, LDL=FALSE, Imult=1)
    rm(retrms, glmod,family)
    if (nAGQ > 0L) {
        function(thetabeta) {
            Lambdat@x[] <<- thfun(thetabeta) # update Lambdat
            LtZt <- Lambdat %*% Zt
            beta[] <<- thetabeta[betaind]
            offb <- offset + X %*% beta
            updatemu <- function(uu) {
                eta[] <- offb + as.vector(crossprod(LtZt, uu))
                mu[] <- linkinv(eta)
                sum(sqDevResid(y, mu, weights)) + sum(uu^2)
            }
            u[] <- 0
            olducden <- updatemu(u)
                                        # initialization
            cvgd <- FALSE
            for(i in 1:npirls){         # PIRLS
                                        # update w and muEta
                Whalf <- Diagonal(x=sqrt(weights/variance(mu)))
                LtZtMWhalf <- LtZt %*% (Diagonal(x=muEta(eta)) %*% Whalf)
                L <- update(L, LtZtMWhalf, 1)
                delu <- as.vector(solve(L, LtZtMWhalf %*%(Whalf %*% (y - mu)) - u))
                if (verbose > 0L) {
                    cat(sprintf("inc: %12.4g", delu[1]))
                    nprint <- min(5,length(delu))
                    for (j in 2:nprint) cat(sprintf(" %12.4g", delu[j]))
                    cat("\n")
                }
                ucden <- updatemu(u + delu)
                if (verbose > 1L) {
                    cat(sprintf("%6.4f: %10.3f\n", 1, ucden))
                }
                
                if(abs((olducden - ucden) / ucden) < tol){
                    cvgd <- TRUE
                    break
                }
                                        # step-halving
                if(ucden > olducden){
                    for(j in 1:10){
                        ucden <- updatemu(u + (delu <- delu/2))
                        if (verbose > 1L) {
                            cat(sprintf("%6.4f: %10.3f\n", 1/2^j, ucden))
                        }
                        if(ucden < olducden) break
                    }
                    if(ucden > olducden) stop("Step-halving failed")
                }
                olducden <- ucden
                u[] <<- u + delu
            }
            if(!cvgd) stop("PIRLS failed to converge")

            ldL2 <- 2*determinant(L, logarithm = TRUE)$modulus
            attributes(ldL2) <- NULL
            ## create the Laplace approx to -2log(L)
            Lm2ll <- aic(y,rep.int(1,n),mu,weights,NULL) + sum(u^2) +
                ldL2 + (q/2)*log(2*pi)
            if (verbose > 0L) {
                cat(sprintf("%10.3f: %12.4g", Lm2ll, thetabeta[1]))
                for (j in 2:length(thetabeta)) cat(sprintf(" %12.4g", thetabeta[j]))
                cat("\n")
            }
            Lm2ll
        }
    } else stop("code for nAGQ == 0 needs to be added")
}
