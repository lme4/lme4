##' Create a univariate Gauss-Hermite quadrature rule
##'
##' This version of Gauss-Hermite quadrature provides the node
##' positions and weights for a scalar integral of a function
##' multiplied by the standard normal density.
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
    ## library("fastGHQuad")
    ## rescale <- function(x, scale.weights=TRUE, scale.roots=TRUE) {
    ##     within(x, {
    ##         if (scale.weights) w <- w/sum(w)
    ##         if (scale.roots) x <- x*sqrt(2)
    ##     })
    ## }
    ## fgqRules <- lapply(1:100, function(n) setNames(rescale(gaussHermiteData(n)), c("z", "w")))
    ## ## However, this shows small differences !!
    ## all.equal(lme4:::fgq_rules, fgqRules, tol=0)

    fr <- as.data.frame(fgq_rules[[ord]])

    rownames(fr) <- NULL
    fr$ldnorm <- dnorm(fr$z, log=TRUE)
    if (asMatrix) as.matrix(fr) else fr
}
