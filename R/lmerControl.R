namedList <- function(...) {
    L <- list(...)
    snm <- sapply(substitute(list(...)),deparse)[-1]
    if (is.null(nm <- names(L))) nm <- snm
    if (any(nonames <- nm=="")) nm[nonames] <- snm[nonames]
    setNames(L,nm)
}
## TESTING:
## a <- b <- c <- 1
## namedList(a,b,c)
## namedList(a,b,d=c)
## namedList(e=a,f=b,d=c)

##' @title Control of mixed model fitting
##' @param optimizer character - name of optimizing function(s).
##'     A character vector or list of functions: length 1 for \code{lmer}
##'    or \code{glmer}, possibly length 2 for \code{glmer}). The built-in optimizers are
##'    \code{\link{Nelder_Mead}} and \code{\link[minqa]{bobyqa}} (from
##'    the \pkg{minqa} package). Any minimizing function that allows
##'    box constraints can be used
##'    provided that it (1) takes input parameters \code{fn} (function
##'    to be optimized), \code{par} (starting parameter values),
##'    \code{lower} (lower bounds) and \code{control} (control parameters,
##'    passed through from the \code{control} argument) and (2)
##'    returns a list with (at least) elements \code{par} (best-fit
##'    parameters), \code{fval} (best-fit function value), \code{conv}
##'    (convergence code) and (optionally) \code{message} (informational
##'    message, or explanation of convergence failure).
##'    Special provisions are made for \code{\link{bobyqa}},
##'    \code{\link{Nelder_Mead}}, and optimizers wrapped in
##'    the \pkg{optimx} package; to use \pkg{optimx} optimizers
##'    (including \code{L-BFGS-B} from base \code{\link{optim}} and
##'    \code{\link{nlminb}}), pass the \code{method} argument to \code{optim}
##'    in the \code{control} argument.
##'    
##'    For \code{glmer}, if \code{length(optimizer)==2}, the first element will be used
##'    for the preliminary (random effects parameters only) optimization, while
##'    the second will be used for the final (random effects plus
##'    fixed effect parameters) phase. 
##' @param sparseX logical - should a sparse model matrix be used for the
##'    fixed-effects terms?  Defaults to \code{FALSE}. Currently inactive.
##' @param restart_edge logical - should the optimizer attempt a restart when it finds a solution at the boundary (i.e. zero random-effect variances or perfect +/-1 correlations)?
##' @param check.numlev.gtr.5 character - rules for checking whether all random effects have >= 5 levels. "ignore": skip the test. "warn": warn if test fails. "stop": throw an error if test fails.
##' @param check.rankZ.gtr.obs character - rules for checking whether the rank of the random effects design matrix Z is greater than the number of observations, indicating possible overfitting.  As for \code{check.numlevel.gtr.5}, with the addition of "warnSmall" and "stopSmall", which run the test only if the dimensions of \code{Z} are <1e6.
##' @param \dots additional arguments to be passed to the nonlinear optimizer (see \code{\link{NelderMead}},
##'    \code{\link[minqa]{bobyqa}})
##' @export
lmerControl <- function(optimizer="Nelder_Mead",
                        restart_edge=TRUE,
                        sparseX=FALSE,
                        check.rankZ.gtr.obs="stopSmall",
                        check.numlev.gtr.5="warning",
                        ...) {
    ## FIXME: is there a better idiom?  match.call() ?
    ## FIXME: check list(...) against formals(get(optimizer)) ?
    namedList(optimizer,
              restart_edge,
              checkControl=
                 namedList(check.rankZ.gtr.obs,
                        check.numlev.gtr.5),
              optControl=list(...))
}

##' @rdname lmerControl
##' @param tolPwrss numeric scalar - the tolerance for declaring convergence in
##'    the penalized iteratively weighted residual sum-of-squares step.
##'    Defaults to 1e-7.
##' @param compDev logical scalar - should compiled code be used for the
##'    deviance evaluation during the optimization of the parameter
##'    estimates?  Defaults to \code{TRUE}.
##' @export
glmerControl <- function(optimizer=c("bobyqa","Nelder_Mead"),
                         restart_edge=TRUE,
                         sparseX=FALSE,
                         check.rankZ.gtr.obs="stopSmall",
                         tolPwrss = 1e-7,
                         compDev = TRUE,
                         ...) {
    if (length(optimizer)==1) {
        optimizer <- replicate(2,optimizer)
    }
    namedList(optimizer,
              restart_edge,
              tolPwrss,
              compDev,
              checkControl=namedList(check.rankZ.gtr.obs),
              optControl=list(...))
}


##' @rdname lmerControl
##' @export
nlmerControl <- function(optimizer="Nelder_Mead",
                         tolPwrss = 1e-10,
                         ...) {
    if (length(optimizer)==1) {
        optimizer <- replicate(2,optimizer)
    }
    namedList(optimizer,
              tolPwrss,
              optControl=list(...))
}


                        
                        
