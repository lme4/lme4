\name{lmerControl}
\title{Control of Mixed Model Fitting}
\alias{glmerControl}
\alias{lmerControl}
\alias{nlmerControl}
\alias{convergence}
\alias{.makeCC}
\description{
  Construct control structures for mixed model fitting.  All arguments
  have defaults, and can be grouped into
  \itemize{
    \item general control parameters, most importantly \code{optimizer},
    further \code{restart_edge}, etc;
    \item model- or data-checking specifications, in short
    \dQuote{checking options}, such as \code{check.nobs.vs.rankZ}, or
    \code{check.rankX} (currently not for \code{nlmerControl});
    \item all the parameters to be passed to the optimizer, e.g.,
    maximal number of iterations, passed via the \code{optCtrl} list
    argument.
  }
}
\usage{
lmerControl(optimizer = "bobyqa",% was "Nelder_Mead" till Dec.2013
    restart_edge = TRUE,
    boundary.tol = 1e-5,
    calc.derivs=TRUE,
    use.last.params=FALSE,
    sparseX = FALSE,
    ## input checking options
    check.nobs.vs.rankZ = "ignore", 
    check.nobs.vs.nlev = "stop",
    check.nlev.gtreq.5 = "ignore",
    check.nlev.gtr.1 = "stop",
    check.nobs.vs.nRE="stop",
    check.rankX = c("message+drop.cols", "silent.drop.cols", "warn+drop.cols",
                    "stop.deficient", "ignore"),
    check.scaleX = c("warning","stop","silent.rescale",
                              "message+rescale","warn+rescale","ignore"),
    check.formula.LHS = "stop",
    ## convergence checking options
    check.conv.grad     = .makeCC("warning", tol = 2e-3, relTol = NULL),
    check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4),
    check.conv.hess     = .makeCC(action = "warning", tol = 1e-6),
    ## optimizer args
    optCtrl = list())

glmerControl(optimizer = c("bobyqa", "Nelder_Mead"),
    restart_edge = FALSE,
    boundary.tol = 1e-5,
    calc.derivs=TRUE,
    use.last.params=FALSE,
    sparseX = FALSE,
    tolPwrss=1e-7,
    compDev=TRUE,
    ## input checking options
    check.nobs.vs.rankZ = "ignore",
    check.nobs.vs.nlev = "stop",
    check.nlev.gtreq.5 = "ignore",
    check.nlev.gtr.1 = "stop",
    check.nobs.vs.nRE="stop",
    check.rankX = c("message+drop.cols", "silent.drop.cols", "warn+drop.cols",
                    "stop.deficient", "ignore"),
    check.scaleX  = "warning",
    check.formula.LHS = "stop",
    ## convergence checking options
    check.conv.grad     = .makeCC("warning", tol = 1e-3, relTol = NULL),
    check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4),
    check.conv.hess     = .makeCC(action = "warning", tol = 1e-6),
    ## optimizer args
    optCtrl = list())

nlmerControl(optimizer = "Nelder_Mead", tolPwrss = 1e-10,
             optCtrl = list())

.makeCC(action, tol, relTol, \dots)
}
\arguments{
  \item{optimizer}{character - name of optimizing
    function(s).  A character vector or list of functions: length 1 for
    \code{lmer} or \code{glmer}, possibly length 2 for \code{glmer}).
    The built-in optimizers are \code{\link{Nelder_Mead}} and
    \code{\link[minqa]{bobyqa}} (from the \pkg{minqa} package). Any
    minimizing function that allows box constraints can be used provided
    that it
    \describe{
      \item{(1)}{takes input parameters \code{fn} (function to be
        optimized), \code{par} (starting parameter values), \code{lower}
        (lower bounds) and \code{control} (control parameters, passed
        through from the \code{control} argument) and}
      \item{(2)}{returns a list with (at least) elements \code{par}
        (best-fit parameters), \code{fval} (best-fit function value),
        \code{conv} (convergence code, equal to zero for
	successful convergence) and (optionally) \code{message}
        (informational message, or explanation of convergence failure).}
    }
    Special provisions are made for \code{\link{bobyqa}},
    \code{\link{Nelder_Mead}}, and optimizers wrapped in the
    \pkg{optimx} package; to use the \pkg{optimx} optimizers (including
    \code{L-BFGS-B} from base \code{\link{optim}} and
    \code{\link{nlminb}}), pass the \code{method} argument to
    \code{optim} in the \code{optCtrl} argument (you may also
    need to load the \code{optimx} package manually using
    \code{\link{library}(optimx)} or \code{\link{require}(optimx)}).

    For \code{glmer}, if \code{length(optimizer)==2}, the first element
    will be used for the preliminary (random effects parameters only)
    optimization, while the second will be used for the final (random
    effects plus fixed effect parameters) phase.  See
    \code{\link{modular}} for more information on these two phases.
  }
  \item{calc.derivs}{logical - compute gradient and Hessian of nonlinear
    optimization solution?}
  \item{use.last.params}{logical - should the last value of the
    parameters evaluated (\code{TRUE}), rather than the value of the
    parameters corresponding to the minimum deviance, be returned?
    This is a "backward bug-compatibility" option; use \code{TRUE}
    only when trying to match previous results.}
  \item{sparseX}{logical - should a sparse model matrix be
    used for the fixed-effects terms?
    Currently inactive.}
  \item{restart_edge}{logical - should the optimizer
    attempt a restart when it finds a solution at the
    boundary (i.e. zero random-effect variances or perfect
    +/-1 correlations)?  (Currently only implemented for
    \code{lmerControl}.)}
  \item{boundary.tol}{numeric - within what distance of
    a boundary should the boundary be checked for a better fit?
    (Set to zero to disable boundary checking.)}
  \item{tolPwrss}{numeric scalar - the tolerance for declaring
    convergence in the penalized iteratively weighted residual
    sum-of-squares step.}
  \item{compDev}{logical scalar - should compiled code be
    used for the deviance evaluation during the optimization
    of the parameter estimates?}
  \item{check.nlev.gtreq.5}{character - rules for
    checking whether all random effects have >= 5 levels.
    See \code{action}.}
  \item{check.nlev.gtr.1}{character - rules for checking
    whether all random effects have > 1 level. See \code{action}.}
  \item{check.nobs.vs.rankZ}{character - rules for
    checking whether the number of observations is greater
    than (or greater than or equal to) the rank of the random
    effects design matrix (Z), usually necessary for
    identifiable variances.  As for \code{action}, with
    the addition of \code{"warningSmall"} and \code{"stopSmall"}, which run
    the test only if the dimensions of \code{Z} are < 1e6.
    \code{nobs > rank(Z)} will be tested for LMMs and GLMMs with
    estimated scale parameters; \code{nobs >= rank(Z)} will be tested
    for GLMMs with fixed scale parameter.
    The rank test is done using the
    \code{method="qr"} option of the \code{\link[Matrix]{rankMatrix}}
    function.
  }
  \item{check.nobs.vs.nlev}{ character - rules for checking whether the
    number of observations is less than (or less than or equal to) the
    number of levels of every grouping factor, usually necessary for
    identifiable variances.  As for \code{action}.
    \code{nobs<nlevels} will be tested for LMMs and GLMMs with estimated
    scale parameters; \code{nobs<=nlevels} will be tested for GLMMs with
    fixed scale parameter.}
  \item{check.nobs.vs.nRE}{character - rules for
    checking whether the number of observations is greater
    than (or greater than or equal to) the number of random-effects
    levels for each term, usually necessary for identifiable variances.
    As for \code{check.nobs.vs.nlev}.}
  \item{check.conv.grad}{rules for checking the gradient of the deviance
    function for convergence.  A list as returned
    by \code{.makeCC}, or a character string with only the action.}
  \item{check.conv.singular}{rules for checking for a singular fit,
    i.e. one where some parameters are on the boundary of the feasible
    space (for example, random effects variances equal to 0 or
    correlations between random effects equal to +/- 1.0);
    as for \code{check.conv.grad} above.}
  \item{check.conv.hess}{rules for checking the Hessian of the deviance
    function for convergence.; as for \code{check.conv.grad}
    above.}
  \item{check.rankX}{character - specifying if \code{\link[Matrix]{rankMatrix}(X)}
    should be compared with \code{ncol(X)} and if columns from the design
    matrix should possibly be dropped to ensure that it has full rank.
    Sometimes needed to make the model identifiable.  The options can be
    abbreviated; the three \code{"*.drop.cols"} options all do drop
    columns, \code{"stop.deficient"} gives an error when the rank is
    smaller than the number of columns where \code{"ignore"} does no
    rank computation, and will typically lead to less easily
    understandable errors, later.}
  \item{check.scaleX}{character - check for problematic scaling of
    columns of fixed-effect model matrix, e.g. parameters measured on
    very different scales.}
  \item{check.formula.LHS}{ check whether specified formula has
    a left-hand side. Primarily for internal use within
    \code{simulate.merMod};
    \emph{use at your own risk} as it may allow the generation
    of unstable \code{merMod} objects}

  \item{optCtrl}{a \code{\link{list}} of additional arguments to be
    passed to the nonlinear optimizer (see \code{\link{Nelder_Mead}},
    \code{\link[minqa]{bobyqa}}).  In particular, both
    \code{Nelder_Mead} and \code{bobyqa} use \code{maxfun} to
    specify the maximum number of function evaluations they
    will try before giving up - in contrast to
    \code{\link{optim}} and \code{optimx}-wrapped optimizers,
    which use \code{maxit}.}
  \item{action}{character - generic choices for the severity level
    of any test. "ignore": skip the test. "warning": warn if test fails.
    "stop": throw an error if test fails.}
  \item{tol}{numeric - tolerance for check }
  \item{relTol}{numeric - tolerance for checking relative variation}
  \item{\dots}{other elements to include in check specification}
}
\value{
  The \code{*Control} functions return a list (inheriting from class
  \code{"merControl"}) containing
  \enumerate{
    \item general control parameters, such as \code{optimizer}, \code{restart_edge};
    \item (currently not for \code{nlmerControl}:)
      \code{"checkControl"}, a \code{\link{list}} of data-checking
      specifications, e.g., \code{check.nobs.vs.rankZ};
    \item parameters to be passed to the optimizer, i.e., the \code{optCtrl}
    list, which may contain \code{maxiter}.
  }

  \code{.makeCC} returns a list containing the check specification
  (action, tolerance, and optionally relative tolerance).
  
}
\details{
  Note that (only!) the pre-fitting \dQuote{checking options}
  (i.e., all those starting with \code{"check."} but \emph{not}
  including the convergence checks (\code{"check.conv.*"}) or
  rank-checking (\code{"check.rank*"}) options)
  may also be set globally via \code{\link{options}}.
  In that case, \code{(g)lmerControl} will use them rather than the
  default values, but will \emph{not} override values that are passed as
  explicit arguments.

  For example, \code{options(lmerControl=list(check.nobs.vs.rankZ = "ignore"))}
  will suppress warnings that the number of observations is less than
  the rank of the random effects model matrix \code{Z}.
}
\section{Troubleshooting}{
  Some additional information about specific warnings and what to do
  about them:
  \itemize{
    \item \code{Model failed to converge with max|grad| ...} The scaled
    gradient at the fitted (RE)ML estimates
    is worryingly large. Try switching optimizers to see if you
    get consistent results, or
    compute values of the deviance in the neighbourhood of the estimated
    parameters to double-check that \code{lme4} has really found a
    local optimum.
    \item \code{Model is nearly unidentifiable ... Rescale variables?}
    The Hessian (curvature matrix) at the maximum likelihood or REML
    estimates has a large eigenvalue, indicating that the surface is
    nearly flat in some direction. Consider centering and/or scaling
    continuous predictor variables.
  }
}
\examples{
str(lmerControl())
str(glmerControl())
\dontrun{
    ## fit with default Nelder-Mead algorithm ...
    fm0 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
    fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
    ## or with minqa::bobyqa ...
    fm1_bobyqa <- update(fm1,control=lmerControl(optimizer="bobyqa"))
    ## or with the nlminb function used in older (<1.0) versions of lme4;
    ## this will usually replicate older results
    require(optimx)
    fm1_nlminb <- update(fm1,control=lmerControl(optimizer="optimx",
                             optCtrl=list(method="nlminb")))
    ## The other option here is method="L-BFGS-B".
    ## Or we can wrap base::optim():
    optimwrap <- function(fn,par,lower,upper,control=list(),
                          ...) {
        if (is.null(control$method)) stop("must specify method in optCtrl")
        method <- control$method
        control$method <- NULL
        ## "Brent" requires finite upper values (lower bound will always
        ##  be zero in this case)
        if (method=="Brent") upper <- pmin(1e4,upper)
        res <- optim(par=par,fn=fn,lower=lower,upper=upper,
                     control=control,method=method,...)
        with(res,list(par=par,
                      fval=value,
                      feval=counts[1],
                      conv=convergence,
                      message=message))
    }
    fm0_brent <- update(fm0,control=lmerControl(optimizer="optimwrap",
                            optCtrl=list(method="Brent")))
    ## You can also use functions from the nloptr package.
    if (require(nloptr)) {
        defaultControl <- list(algorithm="NLOPT_LN_BOBYQA",
                               xtol_abs=1e-6,ftol_abs=1e-6,maxeval=1e5)
        nloptwrap <- function(fn,par,lower,upper,control=list(),...) {
            for (n in names(defaultControl))
                if (is.null(control[[n]])) control[[n]] <- defaultControl[[n]]
            res <- nloptr(x0=par,eval_f=fn,lb=lower,ub=upper,opts=control,...)
            with(res,list(par=solution,
                          fval=objective,
                          feval=iterations,
                          conv=if (status>0) 0 else status,
                          message=message))
        }
        fm1_nloptr <- update(fm1,control=lmerControl(optimizer="nloptwrap"))
        fm1_nloptr_NM <- update(fm1,control=lmerControl(optimizer="nloptwrap",
                                    optCtrl=list(algorithm="NLOPT_LN_NELDERMEAD")))
    }
    ## other algorithm options include NLOPT_LN_COBYLA, NLOPT_LN_SBPLX
  }
}
