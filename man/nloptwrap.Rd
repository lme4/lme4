\name{nloptwrap}
\alias{nloptwrap}
\title{Wrapper for nlOpt optimizers}
\description{
  Wrapper to allow use of \code{nlOpt} optimizers for
  nonlinear optimization stage

}
\usage{
nloptwrap(fn,par,lower,upper,control=list(),...)
}
\arguments{
  \item{fn}{objective function}
  \item{par}{starting parameter vector}
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
\examples{
environment(nloptwrap)$defaultControl
library(lme4)
fm1 <- lmer(Reaction~Days+(Days|Subject),sleepstudy)
fm1_nloptr <- update(fm1,control=lmerControl(optimizer="nloptwrap"))
fm1_nloptr_NM <- update(fm1,control=lmerControl(optimizer="nloptwrap",
                            optCtrl=list(algorithm="NLOPT_LN_NELDERMEAD")))
## other algorithm options include NLOPT_LN_COBYLA, NLOPT_LN_SBPLX
}
