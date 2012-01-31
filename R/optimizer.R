##' Nelder-Mead optimization of parameters that may be subject to box constraints
##'
##' @title Nelder-Mead optimization
##' @param ff a function of one numeric vector argument returning a numeric scalar
##' @param lower numeric vector of lower bounds - elements may be \code{-Inf}.
##' @param upper numeric vector of upper bounds - elements may be \code{Inf}.
##' @param xst numeric vector of initial step sizes to establish the simplex -
##'     all elements must be non-zero.
##' @param x0 numeric vector of starting values for the parameters.
##' @param xt numeric vector of tolerances on the parameters.
##' @param control a named list of control settings.  Possible settings are
##' \describe{
##'     \item{iprint}{numeric scalar - frequency of printing evaluation information.
##'                   Defaults to 0 indicating no printing.}
##'     \item{maxfun}{numeric scalar - maximum number of function evaluations allowed.}
##'     \item{FtolAbs}{numeric scalar - absolute tolerance on change in function values}
##'     \item{FtolRel}{numeric scalar - relative tolerance on change in function values}
##'     \item{XtolRel}{numeric scalar - relative tolerance on change in parameter values}
##'     \item{MinfMax}{numeric scalar - maximum value of the minimum}
##' }
##' @return a list with 4 components
##' \item{fval}{numeric scalar - the minimum function value achieved}
##' \item{pars}{numeric vector - the value of \code{x} providing the minimum}
##' \item{code}{integer scalar - convergence code}
##' \item{control}{list - the list of control settings after substituting for defaults}
##' @export
Nelder_Mead <- function(ff, x0, xst, xt, lower=rep.int(-Inf, n),
                        upper=rep.int(Inf, n), control=list()) {
    stopifnot(is.function(ff),
              length(formals(ff)) == 1L,
              (n <- length(x0 <- as.numeric(x0))) == length(lower <- as.numeric(lower)),
              length(upper <- as.numeric(upper)) == n,
              length(xst <- as.numeric(xst)) == n,
              all(xst != 0),
              length(xt <- as.numeric(xt)) == n)
    nM <- NelderMead$new(lower=lower, upper=upper, x0=x0, xst=xst, xt=xt)
    cc <- do.call(function(iprint=0L, maxfun=10000L, FtolAbs=1e-5,
                           FtolRel=1e-15, XtolRel=1e-7,
                           MinfMax=.Machine$double.xmin, ...) {
        if (length(list(...))>0) warning("unused control arguments ignored")
        list(iprint=iprint, maxfun=maxfun, FtolAbs=FtolAbs, FtolRel=FtolRel,
             XtolRel=XtolRel, MinfMax=MinfMax)
    }, control)
    nM$setFtolAbs(cc$FtolAbs)
    nM$setFtolRel(cc$FtolRel)
    nM$setIprint(cc$iprint)
    nM$setMaxeval(cc$maxfun)
    nM$setMinfMax(cc$MinfMax)
    while ((nMres <- nM$newf(ff(nM$xeval()))) == 0L) {}
    list(fval=nM$value(), pars=nM$xpos(), code=nMres, control=cc)
}
