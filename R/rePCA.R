##' PCA of random-effects variance-covariance estimates
##'
##' Perform a Principal Components Analysis (PCA) of the random-effects
##' variance-covariance estimates from a fitted mixed-effects model
##' @title PCA of random-effects
##' @param x a merMod object
##' @return a \code{prcomplist} object
##' @author Douglas Bates
##' @export
rePCA <- function(x) UseMethod('rePCA')

#' @export
rePCA.merMod <- function(x) {
    chfs <- getME(x,"Tlist")  # list of lower Cholesky factors
    nms <- names(chfs)
    unms <- unique(nms)
    names(unms) <- unms
    svals <- function(m) {
        vv <- svd(m,nv=0L)
        names(vv) <- c("sdev","rotation")
        vv$center <- FALSE
        vv$scale <- FALSE
        class(vv) <- "prcomp"
        vv
    }
    structure(lapply(unms,function(m) svals(Matrix::bdiag(chfs[which(nms == m)]))),
              class="prcomplist")
}
#' @export
summary.prcomplist <- function(object,...) {
    lapply(object,summary)
}
