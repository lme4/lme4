\docType{methods}
\name{profile-methods}
\alias{profile.merMod}
\alias{profile-methods}
\title{Methods for profile() of [ng]lmer fitted models}
\usage{
  \method{profile}{merMod} (fitted, alphamax = 0.01,
    maxpts = 100, delta = cutoff/8, verbose = 0,
    devtol = 1e-09, startval = "prev",
    optimizer = "bobyqa", ...)
}
\arguments{
  \item{fitted}{a fitted model, e.g., the result of
  \code{\link{lmer}(..)}.}

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

  \item{\dots}{potential further arguments for
  \code{profile} methods.}
}
\description{
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
system.time( tpr <- profile(fm01ML) )
(confint(tpr) -> CIpr)
print(xyplot(tpr))
}
\seealso{
  For (more expensive) alternative confidence intervals:
  \code{\link{bootMer}}.
}
\keyword{methods}

