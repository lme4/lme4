## FIX ME: get rid of unnecessary/leftover stuff;
## right now we are only using deriv12??

#' @title Compute gradient of a function.
#'
#' @description
#' This function computes the gradient (first derivative) of a function
#' using a central finite difference method.
#'
#' @details
#' The function has to return a single scalar numerical value and take
#' as argument a numeric vector of the same length as \code{x}. The function also
#' has to be evaluable in a neighborhood around \code{x}, i.e.~at
#' \code{x - delta} and at \code{x + delta} for all elements in \code{x}.
#'
#' @param fun function for which to compute the gradient.
#' @param x numeric vector of values at which to compute the
#' gradient.
#' @param delta the amount to subtract and add to \code{x} when computing the
#' central difference.
#' @param ... additional arguments to \code{fun}.
#'
#' @family derivatives
#' @author Rune Haubo Bojesen Christensen
#' @keywords utilities
#' @export
#' @return a numeric vector of gradients of the same length as
#' \code{x}.
#' @example inst/examples/examples-gradient.R
## gradient <- function(fun, x, delta=1e-4, ...) {
##     sapply(seq_along(x), function(i) {
##         xadd <- xsub <- x
##         xadd[i] <- x[i] + delta
##         xsub[i] <- x[i] - delta
##         (fun(xadd, ...) - fun(xsub, ...)) / (2 * delta)
##     })
## }

gradient <- function(fun, x, delta=1e-4, lower=rep(NA,length(x)),
                      upper=rep(NA,length(x)), ...) {
    xadd <- xsub <- x
    sapply(seq_along(x), function(i) {
        ldelta <- udelta <- delta
        xadd[i] <- x[i] + udelta
        if (!is.na(uu <- upper[i]) && xadd[i]>uu) {
            udelta <- uu-x[i]
            xadd[i] <- uu
        }
        xsub[i] <- x[i] + ldelta
        if (!is.na(ll <- lower[i]) && xsub[i]<ll) {
            ldelta <- x[i]-ll
            xsub[i] <- ll
        }
        (fun(xadd, ...) - fun(xsub, ...)) / (ldelta + udelta)
    })
}

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
##' @example inst/examples/examples-deriv12.R
deriv12 <- function(fun, x, delta=1e-4, fx=NULL,
                    lower=rep(NA,length(x)), upper=rep(NA,length(x)), ...) {
### Compute gradient and Hessian at the same time (to save computing
### time)
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
        H[j, j] <- fadd/udelta[j]^2 - 2 * fx/(udelta[j]*ldelta[j]) +
            fsub/ldelta[j]^2
        g[j] <- (fadd - fsub) / Delta[j]
        ## Off diagonal elements:
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
    list(gradient = g, Hessian = H)
}

##' This function computes the 2nd derivative aka Hessian of a function
##' with vector-valued input and scalar output using a central finite
##' difference method.
##'
##' The function has to return a single scalar numerical value and take
##' as argument a numeric vector of the same length as \code{x}. The function also
##' has to be evaluable in a neighborhood around \code{x}, i.e. at
##' \code{x - delta} and at \code{x + delta} for all elements in \code{x}.
##'
##' @title Compute Hessian with finite difference
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
##' @return the second derivative matrix.
##' @author Rune Haubo Bojesen Christensen
##' @example inst/examples/examples-hessian.R
hessian <- function(fun, x, delta=1e-4, fx=NULL, ...) {
    nx <- length(x)
    fx <- if(!is.null(fx)) fx else fun(x, ...)
    H <- array(NA, dim=c(nx, nx))
    for(j in 1:nx) {
        ## Diagonal elements:
        xadd <- xsub <- x
        xadd[j] <- x[j] + delta
        xsub[j] <- x[j] - delta
        H[j, j] <- (fun(xadd, ...) - 2 * fx +
                    fun(xsub, ...)) / delta^2
        ## Upper triangular (off diagonal) elements:
        for(i in 1:nx) {
            if(i >= j) break
            xaa <- xas <- xsa <- xss <- x
            xaa[c(i, j)] <- x[c(i, j)] + c(delta, delta)
            xas[c(i, j)] <- x[c(i, j)] + c(delta, -delta)
            xsa[c(i, j)] <- x[c(i, j)] + c(-delta, delta)
            xss[c(i, j)] <- x[c(i, j)] - c(delta, delta)
            H[j, i] <- H[i, j] <-
                (fun(xaa, ...) - fun(xas, ...) -
                 fun(xsa, ...) + fun(xss, ...)) /
                     (4 * delta^2)
        }
    }
    H
}

mygrad <-
    function(fun, x, delta = 1e-4,
             method = c("central", "forward", "backward"), ...)
{
    method <- match.arg(method)
    nx <- length(x)
    if(method %in% c("central", "forward")) {
        Xadd <- matrix(rep(x, nx), nrow=nx, byrow=TRUE) + diag(delta, nx)
        fadd <- apply(Xadd, 1, fun, ...)
    }
    if(method %in% c("central", "backward")) {
        Xsub <- matrix(rep(x, nx), nrow=nx, byrow=TRUE) - diag(delta, nx)
        fsub <- apply(Xsub, 1, fun, ...) ## eval.parent perhaps?
    }
    res <- switch(method,
           "forward" = (fadd - fun(x, ...)) / delta,
           "backward" = (fun(x, ...) - fsub) / delta,
           "central" = (fadd - fsub) / (2 * delta)
           )
    res
}

grad.ctr3 <- function(fun, x, delta=1e-4, ...) {
    nx <- length(x)
    Xadd <- matrix(rep(x, nx), nrow=nx, byrow=TRUE) + diag(delta, nx)
    Xsub <- matrix(rep(x, nx), nrow=nx, byrow=TRUE) - diag(delta, nx)
    fadd <- apply(Xadd, 1, fun, ...)
    fsub <- apply(Xsub, 1, fun, ...) ## eval.parent perhaps?
    (fadd - fsub) / (2 * delta)
}

grad.ctr2 <- function(fun, x, delta=1e-4, ...) {
    ans <- x
    for(i in seq_along(x)) {
        xadd <- xsub <- x
        xadd[i] <- x[i] + delta
        xsub[i] <- x[i] - delta
        ans[i] <- (fun(xadd, ...) - fun(xsub, ...)) / (2 * delta)
    }
    ans
}

grad.ctr <- function(fun, x, delta=1e-4, ...) {
    sapply(seq_along(x), function(i) {
        xadd <- xsub <- x
        xadd[i] <- x[i] + delta
        xsub[i] <- x[i] - delta
        (fun(xadd, ...) - fun(xsub, ...)) / (2 * delta)
    })
}

grad.ctr4 <- function(fun, x, delta=1e-4, ...) {
### - checking finiteness of x and fun-values
### - taking care to avoid floating point errors
### - not using h=x*delta rather than h=delta (important for small or
###   large x?)
    if(!all(is.finite(x)))
        stop("Cannot compute gradient: non-finite argument")
    ans <- x ## return values
    for(i in seq_along(x)) {
        xadd <- xsub <- x ## reset fun arguments
        xadd[i] <- x[i] + delta
        xsub[i] <- x[i] - delta
        ans[i] <- (fun(xadd, ...) - fun(xsub, ...)) / (xadd[i] - xsub[i])
### NOTE: xadd[i] - xsub[i] != 2*delta with floating point arithmetic.
    }
    if(!all(is.finite(ans))) {
        warning("cannot compute gradient: non-finite function values occured")
        ans[!is.finite(ans)] <- Inf
    }
    ans
}
