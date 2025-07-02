
##' This function computes the 1st and 2nd derivative (gradient and
##' Hessian) of a function
##' with vector-valued input and scalar output using a central finite
##' difference method.
##'
##' The function has to return a single scalar numerical value and take
##' as argument a numeric vector of the same length as \code{x}. The function also
##' has to be evaluable in a neighborhood around \code{x}, i.e. at
##' \code{x - delta} and at \code{x + delta} for all elements in \code{x}.
##'
##' @title Compute 1st and 2nd derivative
##' @param fun function for which to compute the gradient.
##' @param x numeric vector of values at which to compute the
##' gradient.
##' @param delta the amount to subtract and add to \code{x} when computing the
##' central difference.
##' @param fx optional value of \code{fun(x, ...)}.
##' @param ... additional arguments to \code{fun}.
##'
##' @keywords utilities
##' @export
##' @family derivatives
##'
##' @return a list with components
##' \item{gradient}{the first derivative vector}
##' \item{Hessian}{the second derivative matrix}
##' @author Rune Haubo Bojesen Christensen
##' 
deriv12 <- function(fun, x, delta=1e-4, fx=NULL,
                    lower=rep(NA,length(x)), upper=rep(NA,length(x)), ...) {
### Compute gradient and Hessian simultaneously (to save computing time)
    nx <- length(x)
    if(is.null(fx)) fx <- fun(x, ...)
    stopifnot(length(fx) == 1, nx >= 1)
    H <- array(NA_real_, dim=c(nx, nx))
    g <- numeric(nx)
    xadd <- x + delta
    hasUB <- !missing(upper) && any(active <- !is.na(upper))
    if(hasUB) { # only then may have non-NA
        if((hasUB <- any(active <- active & xadd > upper))) {
            udelta <- ifelse(active, upper-x, delta)
            xadd[active] <- upper[active]
        }
    }
    xsub <- x - delta
    hasLB <- !missing(lower) && any(active <- !is.na(lower))
    if(hasLB) { # only then may have non-NA
        if((hasLB <- any(active <- active & xsub < lower))) {
            ldelta <- ifelse(active, x-lower, delta)
            xsub[active] <- lower[active]
        }
    }
    ## substitute elements of 'mod' vectors into position(s) 'pos'
    ## in base
    spos <- function(base,mod,pos) {
        if (is.list(mod)) {
            for (i in seq_along(mod)) {
                base <- spos(base,mod[[i]],pos[i])
            }
            base
        } else {
            base[pos] <- mod[pos]
            base
        }
    }
    for(j in 1:nx) {
        ## Diagonal elements:
        fadd <- fun(spos(x,xadd,j), ...)
        fsub <- fun(spos(x,xsub,j), ...)
        udj <- if(hasUB) udelta[j] else delta
        ldj <- if(hasLB) ldelta[j] else delta
        H[j, j] <- fadd/udj^2 - 2 * fx/(udj*ldj) + fsub/ldj^2
        g[j] <- (fadd - fsub) / (udj+ldj) #  "total delta" (udelta+ldelta)[j]
        ## Off diagonal elements:
        for(i in seq_len(j - 1L)) { #  1 <= i < j
            ## Compute upper triangular elements (and mirror to lower tri.):
            ud.i <- if(hasUB) udelta[i] else delta
            ld.i <- if(hasLB) ldelta[i] else delta
            xaa <- spos(x,list(xadd,xadd),c(i,j))
            xas <- spos(x,list(xadd,xsub),c(i,j))
            xsa <- spos(x,list(xsub,xadd),c(i,j))
            xss <- spos(x,list(xsub,xsub),c(i,j))
            H[i, j] <- H[j, i] <-
                fun(xaa, ...)/(ud.i+udj)^2 -
                    fun(xas, ...)/(ud.i+ldj)^2 -
                        fun(xsa, ...)/(ld.i+udj)^2 +
                            fun(xss, ...)/(ld.i+ldj)^2
        }
    }
    list(gradient = g, Hessian = H)
}
