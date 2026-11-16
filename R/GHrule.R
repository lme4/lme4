## Create a univariate Gauss-Hermite quadrature rule
##
## This version of Gauss-Hermite quadrature provides the node
## positions and weights for a scalar integral of a function
## multiplied by the standard normal density.
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
    ## Repositioning them similar to how it's done in MixedModels.jl,
    ## by averaging the weights in forward and reverse order and the locations 
    ## in forward order and negative locations in reverse order.
    ## See: https://github.com/lme4/lme4/issues/968#issuecomment-4467149014
    fr$z <- (fr$z - rev(fr$z)) / 2
    fr$w <- (fr$w + rev(fr$w)) / 2
    fr$ldnorm <- dnorm(fr$z, log=TRUE)
    if (asMatrix) as.matrix(fr) else fr
}

## Create a k-dimensional quadrature rule for adaptive Gauss-Hermite
## quadrature in glmer()
##
## Returns a matrix with one node per row: columns 1..k are the node
## coordinates (for a k-dimensional standard normal), column k+1 the
## weight (weights sum to 1) and column k+2 the log standard-normal
## density at the node, summed over dimensions -- the layout expected
## by glmerAGQ() in ../src/external.cpp.  For k == 1 and
## type = "product" this is identical to GHrule(ord).
##
## type = "product": tensor product of the 1-D rule GHrule(ord), i.e.
##   ord^k nodes.  Nodes whose weight is below prune * (largest weight)
##   are dropped (and the remaining weights rescaled to sum to 1).
## type = "sparse": the Smolyak sparse grid GQdk(k, ord), which has far
##   fewer nodes for k > 1 but some negative weights.
AGQrule <- function(ord, k = 1L, type = c("product", "sparse"), prune = 0) {
    type <- match.arg(type)
    stopifnot(length(k) == 1L, (k <- as.integer(k)) >= 1L,
              length(prune) == 1L, prune >= 0, prune < 1)
    if (ord == 0L) return(matrix(0, nrow = 0L, ncol = k + 2L))
    if (type == "product") {
        gh <- GHrule(ord)
        if (k == 1L && prune == 0) return(gh)
        ii <- as.matrix(expand.grid(rep(list(seq_len(nrow(gh))), k),
                                    KEEP.OUT.ATTRS = FALSE))
        z <- matrix(gh[ii, "z"], ncol = k)
        w <- apply(matrix(gh[ii, "w"], ncol = k), 1L, prod)
        if (prune > 0) {
            keep <- w >= prune * max(w)
            z <- z[keep, , drop = FALSE]
            w <- w[keep] / sum(w[keep])
        }
    } else {
        if (k > length(GQN) || ord > length(GQN[[k]]))
            stop(gettextf("sparse grids are only available for k <= %d and nAGQ <= %d",
                          length(GQN), if (k <= length(GQN)) length(GQN[[k]]) else 0L))
        g <- GQdk(k, ord)
        w <- g[, 1L]
        z <- g[, -1L, drop = FALSE]
    }
    ans <- cbind(z, w, rowSums(dnorm(z, log = TRUE)))
    dimnames(ans) <- list(NULL, c(paste0("z", seq_len(k)), "w", "ldnorm"))
    ans
}

## AGQrule() plus a warning when the number of nodes (each costing one
## evaluation of the linear predictor per deviance evaluation) is large
mkAGQmat <- function(nAGQ, k, type = "product", prune = 0, maxNodes = 1000L) {
    GQmat <- AGQrule(nAGQ, k, type = type, prune = prune)
    if (nrow(GQmat) > maxNodes)
        warning(gettextf(paste("nAGQ = %d for a %d-dimensional random effect uses %d quadrature nodes;",
                               "consider a smaller nAGQ, glmerControl(AGQ_grid = \"sparse\")",
                               "or glmerControl(AGQ_prune = )"),
                         nAGQ, k, nrow(GQmat)), call. = FALSE)
    GQmat
}
