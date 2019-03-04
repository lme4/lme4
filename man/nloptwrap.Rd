\name{nloptwrap}
\alias{nloptwrap}
\alias{nlminbwrap}
\title{Wrappers for additional optimizers}
\description{
  Wrappers to allow use of alternative optimizers,
  from \code{NLopt} library or elsewhere, for
  nonlinear optimization stage
}
\usage{
nloptwrap(par, fn, lower, upper, control=list(),...)
nlminbwrap(par, fn, lower, upper, control=list(),...)
}
\arguments{
  \item{par}{starting parameter vector}
  \item{fn}{objective function}
  \item{lower}{vector of lower bounds}
  \item{upper}{vector of upper bounds}
  \item{control}{list of control parameters}
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
(via the \code{nloptr} R package), and to the \code{\link{nlminb}}
optimizer from base R. (\code{nlminb} is also available via the
\code{optimx} package; this wrapper provides access to \code{nlminb}
without the need to install/link the package, and without the additional
post-fitting checks that are implemented by \code{optimx} (see examples
below).

One important difference between the \code{nloptr}-provided
implementation of BOBYQA and the \code{minqa}-provided version accessible
via \code{optimizer="bobyqa"} is that it provides simpler access to
optimization tolerances. \code{minqa::bobyqa} provides only the \code{rhoend}
parameter ("[t]he smallest value of the trust region radius that is
allowed"), while \code{nloptr} provides
a more standard set of tolerances for relative
or absolute change in the objective function
or the parameter values (\code{ftol_rel}, \code{ftol_abs},
\code{xtol_rel}, \code{xtol_abs}).

}
\examples{
environment(nloptwrap)$defaultControl
library(lme4)
## BOBYQA (default): note "nloptwrap" is now the default optimizer
## for lmer
fm1 <- lmer(Reaction~Days+(Days|Subject),
            sleepstudy,
            control=lmerControl(optimizer="nloptwrap"))
## tighten tolerances
fm1B <- update(fm1,control=lmerControl(optimizer="nloptwrap",
                            optCtrl=list(xtol_abs=1e-8,ftol_abs=1e-8)))
## run for longer (no effect in this case)
fm1C <- update(fm1,control=lmerControl(optimizer="nloptwrap",
                                optCtrl=list(maxeval=10000)))

logLik(fm1B)-logLik(fm1)  ## small difference in log likelihood
logLik(fm1C)-logLik(fm1)  ## no difference in LL
## Nelder-Mead
fm1_nloptr_NM <- update(fm1,control=lmerControl(optimizer="nloptwrap",
                            optCtrl=list(algorithm="NLOPT_LN_NELDERMEAD")))
## other nlOpt algorithm options include NLOPT_LN_COBYLA, NLOPT_LN_SBPLX
fm1_nlminb <- update(fm1,control=lmerControl(optimizer="nlminbwrap"))
if (require(optimx)) {
    fm1_nlminb2 <- update(fm1,control=lmerControl(optimizer="optimx",
                              optCtrl=list(method="nlminb",kkt=FALSE)))
}
}
\author{Gabor Grothendieck (\code{nlminbwrap})}
