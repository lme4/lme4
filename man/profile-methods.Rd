\docType{methods}
\name{profile-methods}
\alias{as.data.frame.thpr}
\alias{profile-methods}
\alias{profile.merMod}
\title{Profile method for merMod objects}
\usage{
  \method{profile}{merMod} (fitted, which = 1:nptot,
    alphamax = 0.01, maxpts = 100, delta = cutoff/8,
    verbose = 0, devtol = 1e-09, maxmult = 10,
    startmethod = "prev", optimizer = "bobyqa", ...)

  \method{as.data.frame}{thpr} (x, ...)
}
\arguments{
  \item{fitted}{a fitted model, e.g., the result of
  \code{\link{lmer}(..)}.}

  \item{which}{(integer) which parameters to profile:
  default is all parameters. The parameters are ordered as
  follows: (1) random effects (theta) parameters; (2)
  residual standard deviation (or scale parameter for GLMMs
  where appropriate); (3) fixed effect parameters.  Random
  effects parameters are ordered as in
  \code{getME(.,"theta")}, i.e. as the lower triangle of a
  matrix with standard deviations on the diagonal and
  correlations off the diagonal.}

  \item{alphamax}{maximum alpha value for likelihood ratio
  confidence regions; used to establish the range of values
  to be profiled}

  \item{maxpts}{maximum number of points (in each
  direction, for each parameter) to evaluate in attempting
  to construct the profile}

  \item{delta}{stepping scale for deciding on next point to
  profile}

  \item{verbose}{level of output from internal
  calculations}

  \item{devtol}{tolerance for fitted deviances less than
  baseline (supposedly minimum) deviance}

  \item{maxmult}{maximum multiplier of the original step
  size allowed, defaults to 10.}

  \item{startmethod}{method for picking starting conditions
  for optimization (STUB)}

  \item{optimizer}{(character or function) optimizer to use
  (see \code{\link{lmer}} for details)}

  \item{\dots}{potential further arguments for
  \code{profile} methods.}

  \item{x}{the result of \code{\link{profile}} (or very
  similar structure)}
}
\description{
  Methods for profile() of [ng]lmer fitted models
}
\details{
  Methods for function \code{\link{profile}} (package
  \pkg{stats}), here for profiling (fitted) mixed effect
  models.
}
\section{Methods}{
  FIXME: These (signatures) will change soon --- document
  \bold{after} change! \describe{ \item{signature(fitted =
  \"merMod\")}{ ...  } }
}
\examples{
fm01ML <- lmer(Yield ~ 1|Batch, Dyestuff, REML = FALSE)
## 0.8s (on a 5600 MIPS 64bit fast(year 2009) desktop "AMD Phenom(tm) II X4 925"):
## This is slower because of the Nelder-Mead optimizer but using bobyqa, the default,
## produces a warning.
system.time( tpr  <- profile(fm01ML, optimizer="Nelder_Mead") )
system.time( tpr2 <- profile(fm01ML, optimizer="bobyqa") )
(confint(tpr) -> CIpr)
stopifnot(all.equal(CIpr, confint(tpr2), tol= 1e-11))
xyplot(tpr)
densityplot(tpr, main="densityplot( profile(lmer(..)) )")
splom(tpr)

tpr2 <- profile(fm01ML, which=1:2, optimizer="Nelder_Mead") ## Batch and residual variance only

## GLMM example (running time ~8 seconds on a modern machine)
\dontrun{
gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
            data = cbpp, family = binomial)
system.time(pr4 <- profile(gm1))  ## ~ 7 seconds
xyplot(pr4,layout=c(5,1),as.table=TRUE)
splom(pr4)\% FIXME! currently fails!
}
##'
}
\seealso{
  For (more expensive) alternative confidence intervals:
  \code{\link{bootMer}}.
}
\keyword{methods}

