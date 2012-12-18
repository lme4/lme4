#' Create an RSC object from the rv and xv components
#' 
#' @param fl the matrix of row indices for the regular sparse column representation of Zt
#' @param xv the non-zero values in ZtXt
#' @param theta an optional variance-component parameter vector
#' @param lower optional lower bounds on the variance-component parameter vector
#' @return an RSC object
#' @examples
#' pred <- createRSC(rv=as.integer(Dyestuff$Batch),
#'                   xv=matrix(1, nrow=2L, ncol=nrow(Dyestuff)))
#' str(pred)
#' @export
createRSC <- function(fl, xv, theta=rep.int(1, k), lower=numeric(k)) {
    stopifnot(is.matrix(xv), is.double(xv))
    ## fl can be specified as 0-based or 1-based indices
    stopifnot((minfl <- min(fl <- as.integer(fl))) %in% 0:1)
    if (minfl == 1L) fl <- fl - 1L
    n <- ncol(xv)
    ## regenerating fl as a matrix allows passing a vector in simple cases
    fl <-matrix(fl, ncol=n) 
    q <- max(fl) + 1L
    k <- nrow(fl)
    kpp <- nrow(xv)
    p <- kpp - k
    qpp <- q + p
    stopifnot(p > 0L) # perhaps >= ? Would p == 0 ever make sense?
    theta <- as.numeric(theta)
    lower <- as.numeric(lower)
    stopifnot(length(theta) == length(lower),
              sum(is.finite(lower)) == k)
    i <- do.call(rbind, c(list(fl), as.list(q:(qpp - 1L))))
    A <- tcrossprod(sparseMatrix(i=as.vector(i),
                                 j=rep(0:(n-1), each=nrow(xv)),
                                 x=as.vector(xv), index1=FALSE)) +
                                     Diagonal(x=rep.int(c(1,0), c(q,p)))
    ## Cholesky results are not saved in this function but are cached as part of the A object
    # need to work out the permutation keeping Z and X parts distinct. For now set perm=FALSE
    Cholesky(A, perm=FALSE, LDL=FALSE) 
    ubeta <- numeric(qpp)
    new("RSC", x=xv, i=i, theta=theta, lower=lower, A=A, ubeta=numeric(qpp))
}

#' Update for the penalized least squares problem
#' 
#' @param pred an RSC object
#' @param resid current residual
#' @return the 'ldL2' value for the updated A
#' @examples
#' pred <- createRSC(rv=as.integer(Dyestuff$Batch),
#'                   xv=matrix(1, nrow=2L, ncol=nrow(Dyestuff)))
#' str(pred)
#' with(Dyestuff, RSCupdate(pred, Yield))
#' str(pred)
#' @export
RSCupdate <- function(pred, resid) .Call(lme4_RSCupdate, pred, resid)


