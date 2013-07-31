##' @rdname NelderMead
##' @title Nelder-Mead optimization of parameters, possibly with box constraints
##' @param fn a function of a single numeric vector argument returning a numeric scalar
##' @param par numeric vector of starting values for the parameters.
##' @param lower numeric vector of lower bounds (elements may be \code{-Inf}).
##' @param upper numeric vector of upper bounds (elements may be \code{Inf}).
##' @param control a named list of control settings.  Possible settings are
##' \describe{
##'     \item{iprint}{numeric scalar - frequency of printing evaluation information.
##'                   Defaults to 0 indicating no printing.}
##'     \item{maxfun}{numeric scalar - maximum number of function evaluations allowed (default:10000).}
##'     \item{FtolAbs}{numeric scalar - absolute tolerance on change in function values (default: 1e-5)}
##'     \item{FtolRel}{numeric scalar - relative tolerance on change in function values (default:1e-15)}
##'     \item{XtolRel}{numeric scalar - relative tolerance on change in parameter values (default: 1e-7)}
##'     \item{MinfMax}{numeric scalar - maximum value of the minimum (default: .Machine$double.xmin)}
##'     \item{xst}{numeric vector of initial step sizes to establish the simplex -
##'     all elements must be non-zero (default: rep(0.02,length(par)))}
##'     \item{xt}{numeric vector of tolerances on the parameters (default: xst*5e-4)}
##'     \item{verbose}{numeric value: 0=no printing, 1=print every 20 evaluations,
##'         2=print every 10 evalutions, 3=print every evaluation.  Sets \sQuote{iprint},
##'         if specified, but does not override it.}
##' }
##'
##' @return a list with 4 components
##' \item{fval}{numeric scalar - the minimum function value achieved}
##' \item{par}{numeric vector - the value of \code{x} providing the minimum}
##' \item{ierr}{integer scalar - error code (see below)}
##' \item{control}{list - the list of control settings after substituting for defaults}
##' @note
##' Return error codes (\code{ierr}):
##' \describe{
##' \item{-4}{\code{nm_evals}: maximum evaluations reached}
##' \item{-3}{\code{nm_forced}: ?}
##' \item{-2}{\code{nm_nofeasible}: cannot generate a feasible simplex}
##' \item{-1}{\code{nm_x0notfeasible}: initial x is not feasible (?)}
##' }
##' @export
Nelder_Mead <- function(fn, par, lower=rep.int(-Inf, n),
                        upper=rep.int(Inf, n), control=list()) {
    n <- length(par)
    if (is.null(xst <- control[["xst"]])) xst <- rep.int(0.02,n)
    if (is.null(xt <- control[["xt"]])) xt <- xst*5e-4

    control[["xst"]] <- control[["xt"]] <- NULL

    ## mapping between simpler 'verbose' setting (0=no printing, 1=20, 2=10, 3=1)
    ##  and internal 'iprint' control (frequency of printing)
    if (is.null(verbose <- control[["verbose"]])) verbose <- 0
    control[["verbose"]] <- NULL
    if (is.null(control[["iprint"]])) {
      control[["iprint"]] <- switch(as.character(min(as.numeric(verbose),3L)),
                                    "0"=0, "1"=20,"2"=10,"3"=1)
    }
    stopifnot(is.function(fn),
              length(formals(fn)) == 1L,
              (n <- length(par <- as.numeric(par))) == length(lower <- as.numeric(lower)),
              length(upper <- as.numeric(upper)) == n,
              length(xst <- as.numeric(xst)) == n,
              all(xst != 0),
              length(xt <- as.numeric(xt)) == n)
    nM <- NelderMead$new(lower=lower, upper=upper, x0=par, xst=xst, xt=xt)
    cc <- do.call(function(iprint=0L, maxfun=10000L, FtolAbs=1e-5,
                           FtolRel=1e-15, XtolRel=1e-7,
                           MinfMax=-.Machine$double.xmax, ...) {
        if (length(list(...))>0) warning("unused control arguments ignored")
        list(iprint=iprint, maxfun=maxfun, FtolAbs=FtolAbs, FtolRel=FtolRel,
             XtolRel=XtolRel, MinfMax=MinfMax)
    }, control)
    nM$setFtolAbs(cc$FtolAbs)
    nM$setFtolRel(cc$FtolRel)
    nM$setIprint(cc$iprint)
    nM$setMaxeval(cc$maxfun)
    nM$setMinfMax(cc$MinfMax)
    while ((nMres <- nM$newf(fn(nM$xeval()))) == 0L) {}

    cmsg <- "reached max evaluations"
    if (nMres==-4) {
      ## map max evals from error to warning
      cmsg <- warning(sprintf("failure to converge in %d evaluations",cc$maxfun))
      nMres <- 4
    }

    msgvec <- c("nm_forced","cannot generate a feasible simplex","initial x is not feasible",
                "active","minf_max","fcvg","xcvg",  ## FIXME: names (see NelderMead_newf in external.cpp)
                cmsg)

    if (nMres<0) stop(msgvec[nMres+4])

    cc <- c(cc,xst=xst,xt=xt)
    list(fval=nM$value(), par=nM$xpos(), convergence=pmin(nMres,0), message=msgvec[nMres+4],
         control=cc)
}
