\name{nloptwrap}
\alias{nloptwrap}
\alias{nlminbwrap}
\title{Wrappers for additional optimizers}
\description{
  Wrappers to allow use of alternative optimizers, from the \command{NLopt}
  library (via \CRANpkg{nloptr}) or elsewhere, for the nonlinear optimization stage.
}
\usage{
nloptwrap (par, fn, lower, upper, control = list(), ...)
nlminbwrap(par, fn, lower, upper, control = list(), ...)
}
\arguments{
  \item{par}{starting parameter vector}
  \item{fn}{objective function}
  \item{lower, upper}{numeric vector of lower and upper bounds.}
  \item{control}{\code{\link{list}} of control parameters, corresponding to
    \code{optCtrl = *}, e.g., in \code{\link{lmerControl}()}.  For
    \code{nloptwrap}, see \code{defaultControl} in \sQuote{Details}.}
  \item{\dots}{additional arguments to be passed to objective function}
}
\value{
  \item{par}{estimated parameters}
  \item{fval}{objective function value at minimum}
  \item{feval}{number of function evaluations}
  \item{conv}{convergence code (0 if no error)}
  \item{message}{convergence message}
}
\details{
Using alternative optimizers is an important trouble-shooting
tool for mixed models.  These wrappers provide convenient access to
the optimizers provided by Steven Johnson's \code{NLopt} library
(via the \CRANpkg{nloptr} R package), and to the \code{\link{nlminb}}
optimizer from base R.  \code{nlminb} is also available via the
\CRANpkg{optimx} package; this wrapper provides access to \code{nlminb()}
without the need to install/link the package, and without the additional
post-fitting checks that are implemented by \CRANpkg{optimx} (see examples
below).

  One important difference between the \CRANpkg{nloptr}-provided
  implementation of BOBYQA and the \CRANpkg{minqa}-provided version
  accessible via \code{optimizer="bobyqa"} is that it provides simpler
  access to optimization tolerances.  \code{\link[minqa]{bobyqa}} provides
  only the \code{rhoend} parameter (\dQuote{[t]he smallest value of the
  trust region radius that is allowed}), while \pkg{nloptr} provides a more
  standard set of tolerances for relative or absolute change in the
  objective function or the parameter values (\code{ftol_rel},
  \code{ftol_abs}, \code{xtol_rel}, \code{xtol_abs}).

  The default (empty) \code{control} list corresponds to the following settings:
  \describe{
    \item{\code{nlminbwrap}:}{\code{control} exactly corresponds to
      \code{\link{nlminb}()}'s defaults, see there.}
    \item{\code{nloptwrap}:}{\code{environment(nloptwrap)$defaultControl}
      contains the defaults, notably \code{algorithm = "NLOPT_LN_BOBYQA"}.

      \code{nloptr::nloptr.print.options()} shows and explains the many
      possible algorithm and options.}
  }
}

\examples{
library(lme4)
ls.str(environment(nloptwrap)) # 'defaultControl' algorithm "NLOPT_LN_BOBYQA"
## Note that  'optimizer =  "nloptwrap"' is now the default for lmer() :
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
## tighten tolerances
fm1B <- update(fm1, control= lmerControl(optCtrl = list(xtol_abs=1e-8, ftol_abs=1e-8)))
## run for longer (no effect in this case)
fm1C <- update(fm1,control = lmerControl(optCtrl = list(maxeval=10000)))

  logLik(fm1B) - logLik(fm1)  ## small difference in log likelihood
c(logLik(fm1C) - logLik(fm1)) ## no difference in LL
## Nelder-Mead
fm1_nloptr_NM <- update(fm1, control=
                  lmerControl(optimizer = "nloptwrap",
                              optCtrl = list(algorithm = "NLOPT_LN_NELDERMEAD")))
## other nlOpt algorithm options include NLOPT_LN_COBYLA, NLOPT_LN_SBPLX, see
if(interactive())% (the package *is* installed w/ 'lme4')
  nloptr::nloptr.print.options()

fm1_nlminb <- update(fm1, control=lmerControl(optimizer = "nlminbwrap"))
if (require(optimx)) { ## the 'optimx'-based nlminb :
  fm1_nlminb2 <- update(fm1, control=
                lmerControl(optimizer = "optimx",
                            optCtrl = list(method="nlminb", kkt=FALSE)))
  cat("Likelihood difference (typically zero):  ",
      c(logLik(fm1_nlminb) - logLik(fm1_nlminb2)), "\n")
}


}
\author{Gabor Grothendieck (\code{nlminbwrap})}% Ben Bolker for nloptwrap
