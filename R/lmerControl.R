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
##'    fixed effect parameters) phase. See \code{\link{modular}} for more information on
##'    these two phases.
##' @param sparseX logical - should a sparse model matrix be used for the
##'    fixed-effects terms?  Defaults to \code{FALSE}. Currently inactive.
##' @param restart_edge logical - should the optimizer attempt a restart when it finds a solution at the boundary (i.e. zero random-effect variances or perfect +/-1 correlations)?
##' @param check.nlev.gtreq.5 character - rules for checking whether all random effects have >= 5 levels. "ignore": skip the test. "warn": warn if test fails. "stop": throw an error if test fails.
##' @param check.nlev.gtr.1 character - rules for checking whether all random effects have > 1 level. As for \code{check.nlevel.gtr.5}.
##' @param check.nobs.vs.rankZ character - rules for checking whether the number of observations is greater than (or greater than or equal to) the rank of the random effects design matrix (Z), usually necessary for identifiable variances.  As for \code{check.nlevel.gtreq.5}, with the addition of "warnSmall" and "stopSmall", which run the test only if the dimensions of \code{Z} are <1e6. \code{nobs>rank(Z)} will be tested for LMMs and GLMMs with estimated scale parameters; \code{nobs>=rank(Z)} will be tested for GLMMs with fixed scale parameter.
## FIXME: add (something like) allow.nobs.eq.nlev (and/or rankZ), and
## add better logic for distinguishing whether the scale parameter is
## being estimated
##' @param check.nobs.vs.nlev character - rules for checking whether the number of observations is less than (or less than or equal to) the number of levels of every grouping factor, usually necessary for identifiable variances.  As for \code{check.nlevel.gtreq.5}: \code{nobs<nlevels} will be tested for LMMs and GLMMs with estimated scale parameters; \code{nobs<=nlevels} will be tested for GLMMs with fixed scale parameter.
##' @param \dots additional arguments to be passed to the nonlinear optimizer (see \code{\link{Nelder_Mead}},
##'    \code{\link[minqa]{bobyqa}}). In particular, both \code{Nelder_Mead} and \code{bobyqa} use \code{maxfun} to specify
##'    the maximum number of function evaluations they will try before giving up - in contrast to \code{\link{optim}} and \code{optimx}-wrapped optimizers, which use \code{maxit}.
##' @return a list (of class \code{merControl}) containing (1) general control parameters (e.g. \code{optimizer}, \code{restart_edge}); (2) a list of data-checking specifications (e.g. \code{check.nobs.vs.rankZ}); (3) parameters to be passed to the optimizer (i.e., the contents of \dots, for example \code{maxiter})
##' @details if options are set via \code{\link{options}}, [gn]lmerControl will use them rather than the default values (but will not override values that are passed as explicit arguments); for example, \code{options(check.nlev.gtreq.5="ignore")} will suppress warnings that there an insufficient random effects levels for reliable estimation.
##' @export
lmerControl <- function(optimizer="Nelder_Mead",
                        restart_edge=TRUE,
                        sparseX=FALSE,
                        check.nobs.vs.rankZ="warningSmall",
                        check.nobs.vs.nlev="stop",
                        check.nlev.gtreq.5="ignore",
                        check.nlev.gtr.1="stop",
                        optCtrl = list())
{
    ## FIXME: is there a better idiom?  match.call() ?
    ## fill in values from options, but **only if not specified explicitly in arguments**
    ##  (ugh ... is there a better way to do this?  mapply() is clunky:
    ##  http://stackoverflow.com/questions/16276667/using-apply-with-assign-in-r
    stopifnot(is.list(optCtrl))
    if (!is.null(lmerOpts <- getOption("lmerControl"))) {
        for (arg in names(lmerOpts)) {
            if (do.call(missing,list(arg))) ## missing from explicit arguments
                assign(arg,lmerOpts[[arg]])
        }
    }
    structure(namedList(optimizer, restart_edge,
                        checkControl =
		   namedList(check.nobs.vs.rankZ,
                             check.nobs.vs.nlev,
			     check.nlev.gtreq.5,
			     check.nlev.gtr.1),
                        optCtrl=optCtrl),
	      class = c("lmerControl", "merControl"))
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
                         restart_edge=FALSE,
                         sparseX=FALSE,
                         check.nobs.vs.rankZ="warningSmall",
                         check.nobs.vs.nlev="stop",
                         check.nlev.gtreq.5="ignore",
                         check.nlev.gtr.1="stop",
                         tolPwrss = 1e-7,
                         compDev = TRUE,
                         optCtrl = list())
{
    ## FIXME: should try to modularize/refactor/combine with lmerControl if possible
    ## but note different defaults
    ##                lmer        glmer
    ## optimizer    Nelder_Mead  c(Nelder_Mead,bobyqa)
    ## tolPwrss     N/A          1e-7
    ## compDev      N/A          TRUE
    ##
    ## (and possible future divergence)
    stopifnot(is.list(optCtrl))
    if (length(optimizer)==1) {
	optimizer <- replicate(2,optimizer) # works evevn when optimizer is function
    }
    if (!is.null(glmerOpts <- getOption("glmerControl"))) {
        for (arg in names(glmerOpts)) {
            if (do.call(missing,list(arg))) ## missing from explicit arguments
                assign(arg,glmerOpts[[arg]])
        }
    }
    structure(namedList(optimizer,
			restart_edge,
			tolPwrss,
			compDev,
			checkControl=
			namedList(check.nobs.vs.rankZ,
                                  check.nobs.vs.nlev,
				  check.nlev.gtreq.5,
				  check.nlev.gtr.1),
                        optCtrl=optCtrl),
	      class = c("glmerControl", "merControl"))
}


##' @rdname lmerControl
##' @export
nlmerControl <- function(optimizer="Nelder_Mead",
                         tolPwrss = 1e-10,
                         optCtrl = list())
{
    stopifnot(is.list(optCtrl))
    if (length(optimizer)==1) {
        optimizer <- replicate(2,optimizer)
    }
    structure(namedList(optimizer,
			tolPwrss,
			optCtrl=optCtrl),
	      class = c("nlmerControl", "merControl"))
}
