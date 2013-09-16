\name{lmerControl}
\title{Control of Mixed Model Fitting}
\alias{glmerControl}
\alias{lmerControl}
\alias{nlmerControl}
\description{
  Construct a control structure for mixed model fitting.
}
\usage{
lmerControl(optimizer = "Nelder_Mead",
    restart_edge = TRUE, sparseX = FALSE,
    check.nobs.vs.rankZ = "warningSmall",
    check.nobs.vs.nlev = "stop",
    check.nlev.gtreq.5 = "ignore",
    check.nlev.gtr.1 = "stop", optCtrl = list())

glmerControl(optimizer = c("bobyqa", "Nelder_Mead"),
    restart_edge = FALSE, sparseX = FALSE,
    check.nobs.vs.rankZ = "warningSmall",
    check.nobs.vs.nlev = "stop",
    check.nlev.gtreq.5 = "ignore",
    check.nlev.gtr.1 = "stop", tolPwrss = 1e-07,
    compDev = TRUE, optCtrl = list())

nlmerControl(optimizer = "Nelder_Mead", tolPwrss = 1e-10,
             optCtrl = list())
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
        \code{conv} (convergence code) and (optionally) \code{message}
        (informational message, or explanation of convergence failure).}
    }
    Special provisions are made for \code{\link{bobyqa}},
    \code{\link{Nelder_Mead}}, and optimizers wrapped in the
    \pkg{optimx} package; to use \pkg{optimx} optimizers (including
    \code{L-BFGS-B} from base \code{\link{optim}} and
    \code{\link{nlminb}}), pass the \code{method} argument to
    \code{optim} in the \code{control} argument.

    For \code{glmer}, if \code{length(optimizer)==2}, the first element
    will be used for the preliminary (random effects parameters only)
    optimization, while the second will be used for the final (random
    effects plus fixed effect parameters) phase.  See
    \code{\link{modular}} for more information on these two phases.
  }
  \item{sparseX}{logical - should a sparse model matrix be
    used for the fixed-effects terms?  Defaults to \code{FALSE}.
    Currently inactive.}
  \item{restart_edge}{logical - should the optimizer
    attempt a restart when it finds a solution at the
    boundary (i.e. zero random-effect variances or perfect
    +/-1 correlations)?}

  \item{check.nlev.gtreq.5}{character - rules for
    checking whether all random effects have >= 5 levels.
    "ignore": skip the test. "warning": warn if test fails.
    "stop": throw an error if test fails.}

  \item{check.nlev.gtr.1}{character - rules for checking
    whether all random effects have > 1 level. As for
    \code{check.nlevel.gtr.5}.}

  \item{check.nobs.vs.rankZ}{character - rules for
    checking whether the number of observations is greater
    than (or greater than or equal to) the rank of the random
    effects design matrix (Z), usually necessary for
    identifiable variances.  As for \code{check.nlevel.gtreq.5}, with
    the addition of \code{"warningSmall"} and \code{"stopSmall"}, which run
    the test only if the dimensions of \code{Z} are < 1e6.
    \code{nobs > rank(Z)} will be tested for LMMs and GLMMs with
    estimated scale parameters; \code{nobs >= rank(Z)} will be tested
    for GLMMs with fixed scale parameter.}

  \item{check.nobs.vs.nlev}{ character - rules for checking whether the
  number of observations is less than (or less than or equal to) the
  number of levels of every grouping factor, usually necessary for
  identifiable variances.  As for \code{check.nlevel.gtreq.5}.
  \code{nobs<nlevels} will be tested for LMMs and GLMMs with estimated
  scale parameters; \code{nobs<=nlevels} will be tested for GLMMs with
  fixed scale parameter.}

  \item{optCtrl}{a \code{\link{list}} of additional arguments to be
    passed to the nonlinear optimizer (see \code{\link{Nelder_Mead}},
    \code{\link[minqa]{bobyqa}}).  In particular, both
    \code{Nelder_Mead} and \code{bobyqa} use \code{maxfun} to
    specify the maximum number of function evaluations they
    will try before giving up - in contrast to
    \code{\link{optim}} and \code{optimx}-wrapped optimizers,
    which use \code{maxit}.}
  \item{tolPwrss}{numeric scalar - the tolerance for declaring
    convergence in the penalized iteratively weighted residual
    sum-of-squares step.  Defaults to 1e-7.}
  \item{compDev}{logical scalar - should compiled code be
    used for the deviance evaluation during the optimization
    of the parameter estimates?  Defaults to \code{TRUE}.}
}
\value{
  a list (of class \code{"merControl"}) containing
  (1) general control parameters (e.g. \code{optimizer},
  \code{restart_edge});
  (2) a list of data-checking specifications
  (e.g. \code{check.nobs.vs.rankZ});
  (3) parameters to be passed to the optimizer, i.e., the \code{optCtrl}
  list, which may contain \code{maxiter}.)
}
\details{
  If checking options are set via \code{\link{options}}, \code{[g]lmerControl}
  will use them rather than the default values (but will \emph{not}
  override values that are passed as explicit arguments); for example,
  \code{options(lmerControl=list(check.nobs.vs.rankZ = "ignore"))} will
  suppress warnings that the number of observations is less than
  the rank of the random effects model matrix \code{Z}.
}
\examples{
str(lmerControl())
}
