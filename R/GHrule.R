##' Create a univariate Gauss-Hermite quadrature rule
##'
##' This version of Gauss-Hermite quadrature provides the node
##' positions and weights for a scalar integral of a function
##' multiplied by the standard normal density.
##' @title Univariate Gauss-Hermite quadrature rule
##' @param ord scalar integer between 1 and 25 - the order, or number of
##'    nodes and weights, in the rule.  When the function being
##'    multiplied by the standard normal density is a polynomial of
##'    order 2k-1 the rule of order k integrates the product exactly.
##' @param asMatrix logical scalar - should the result be returned as
##'    a matrix.  If \code{FALSE} a data frame is returned.  Defaults
##'    to \code{TRUE}.
##' @return a matrix with \code{ord} rows and three columns which are
##'    \code{z} the node positions, \code{w} the weights and
##'    \code{ldnorm}, the logarithm of the normal density evaluated at
##'    the nodes.
##' @examples
##' (r5 <- GHrule(5, asMatrix=FALSE))
##' ## second, fourth, sixth, eighth and tenth central moments of the
##' ## standard Gaussian density
##' with(r5, sapply(seq(2, 10, 2), function(p) sum(w * z^p)))
##' }
##' @export
GHrule <- function (ord, asMatrix=TRUE) {
    stopifnot(length(ord) == 1,
              (ord <- as.integer(ord)) >= 0L,
              ord < 101L)
    if (ord == 0L) {
        if (asMatrix) return(matrix(0, nrow=0L, ncol=3L))
        stop ("combination of ord==0 and asMatrix==TRUE not implemented")
    }
    ## fgq_rules comes from sysdata.rda:
    ## result of
    ##
    ## library("fastGHQuad") ## version 0.2
    ## rescale <- function(x,scale.weights=TRUE,
    ##                     scale.roots=TRUE) {
    ##     x <- within(x,
    ##             {
    ##                 if (scale.weights) w <- w/sum(w)
    ##                 if (scale.roots) x <- x*sqrt(2)
    ##                 })
    ##     return(x)
    ## }
    ## lapply(1:100,rename(rescale(gaussHermiteData(x)),c(x="z")))

    fr <- as.data.frame(fgq_rules[[ord]])

    rownames(fr) <- NULL
    fr$ldnorm <- dnorm(fr$z, log=TRUE)
    if (asMatrix) as.matrix(fr) else fr
}
