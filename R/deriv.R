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
##' @return alist with components
##' \item{gradient}{the first derivative vector}
##' \item{Hessian}{the second derivative matrix}
##' @author Rune Haubo Bojesen Christensen
deriv12 <- function(fun, x, delta=1e-4, fx=NULL,
                    lower=rep(NA,length(x)), upper=rep(NA,length(x)),
                    calc.hess=TRUE, ...) {
    ## Compute gradient and Hessian at the same time (to save computing time)
    nx <- length(x)
    fx <- if(!is.null(fx)) fx else fun(x, ...)
    stopifnot(length(fx) == 1)
    H <- array(NA, dim=c(nx, nx))
    g <- numeric(nx)
    xadd <- x + delta
    ubActive <- !is.na(upper) & xadd>upper
    udelta <- ifelse(ubActive,upper-x,delta)
    xadd[ubActive] <- upper[ubActive]
    xsub <- x - delta
    lbActive <- !is.na(lower) & xadd<lower
    ldelta <- ifelse(lbActive,x-lower,delta)
    xsub[lbActive] <- lower[lbActive]
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
    ## TOTAL delta
    Delta <- ldelta+udelta
    for(j in 1:nx) {
        ## Diagonal elements:
        fadd <- fun(spos(x,xadd,j), ...)
        fsub <- fun(spos(x,xsub,j), ...)
        if (calc.hess) 
            H[j, j] <- fadd/udelta[j]^2 - 2 * fx/(udelta[j]*ldelta[j]) +
                                          fsub/ldelta[j]^2
        g[j] <- (fadd - fsub) / Delta[j]
        ## Off diagonal elements:
        if (calc.hess) {
            for(i in 1:nx) {
                if(i >= j) break
                ## Compute upper triangular elements:
                xaa <- spos(x,list(xadd,xadd),c(i,j))
                xas <- spos(x,list(xadd,xsub),c(i,j))
                xsa <- spos(x,list(xsub,xadd),c(i,j))
                xss <- spos(x,list(xsub,xsub),c(i,j))
                H[i, j] <- H[j, i] <-
                    fun(xaa, ...)/(udelta[i]+udelta[j])^2 -
                    fun(xas, ...)/(udelta[i]+ldelta[j])^2 -
                        fun(xsa, ...)/(ldelta[i]+udelta[j])^2 +
                            fun(xss, ...)/(ldelta[i]+ldelta[j])^2
            }
        }
    }
    list(gradient = g, Hessian = H)
}
