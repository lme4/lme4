\docType{methods}
\name{profile-methods}
\alias{as.data.frame.thpr}
\alias{profile-methods}
\alias{profile.merMod}
\title{Profile method for merMod objects}
\usage{
\method{profile}{merMod}(fitted, which = 1:nptot, alphamax = 0.01,
	maxpts = 100, delta = cutoff/8, verbose = 0, devtol = 1e-09,
        maxmult = 10, startmethod = "prev", optimizer = "bobyqa",
	signames = TRUE, ...)

\method{as.data.frame}{thpr} (x, ...)
}
\arguments{
  \item{fitted}{a fitted model, e.g., the result of \code{\link{lmer}(..)}.}
  \item{which}{integer indicating which parameters to profile:
    default is all parameters.  The parameters are ordered as
    follows:
    \describe{
      \item{(1)}{random effects (theta) parameters;}
      \item{(2)}{residual standard deviation (or scale parameter for GLMMs
        where appropriate);}
      \item{(3)}{fixed effect parameters.  Random effects parameters are
	ordered as in \code{getME(.,"theta")}, i.e., as the lower
	triangle of a matrix with standard deviations on the diagonal
	and correlations off the diagonal.}
    }
    }%% FIXME: allow parameter names.
  \item{alphamax}{maximum alpha value for likelihood ratio confidence
    regions; used to establish the range of values to be profiled.}
  \item{maxpts}{maximum number of points (in each direction, for each
    parameter) to evaluate in attempting to construct the profile.}
  \item{delta}{stepping scale for deciding on next point to profile.}
  \item{verbose}{level of output from internal calculations.}
  \item{devtol}{tolerance for fitted deviances less than
    baseline (supposedly minimum) deviance.}
  \item{maxmult}{maximum multiplier of the original step size allowed,
    defaults to 10.}
  \item{startmethod}{method for picking starting conditions for
    optimization (STUB).}
  \item{optimizer}{(character or function) optimizer to use (see
    \code{\link{lmer}} for details).}
  \item{signames}{logical indicating if abbreviated names of the form
    \code{.sigNN} should be used; otherwise, names are more meaningful
    (but longer) of the form \code{(sd|cor)_(effects)|(group)}.  Note
    that some code for profile transformations (e.g.,
    \code{\link{varianceProf}}) depends on \code{signames==TRUE}.}
  \item{\dots}{potential further arguments for \code{profile} methods.}
  \item{x}{the result of \code{\link{profile}} (or very similar structure)}
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
  \describe{ \item{signature(fitted = \"merMod\")}{ ...  }
  }
}
\seealso{
  For (more expensive) alternative confidence intervals:
  \code{\link{bootMer}}.
}
\examples{
fm01ML <- lmer(Yield ~ 1|Batch, Dyestuff, REML = FALSE)
system.time( tpr  <- profile(fm01ML, optimizer="Nelder_Mead") )
## ~2.6s (on a 2010 Macbook Pro)
system.time( tpr  <- profile(fm01ML))
## ~1s, + possible warning about bobyqa convergence
(confint(tpr) -> CIpr)
\donttest{% too much precision (etc). but just FYI:
stopifnot(all.equal(CIpr,
  array(c(12.1985292, 38.2299848, 1486.4515,
          84.0630513, 67.6576964, 1568.54849), dim = 3:2,
        dimnames = list(c(".sig01", ".sigma", "(Intercept)"),
                        c("2.5 \%", "97.5 \%"))),
                    tol= 1e-07))# 1.37e-9 {64b}
}
xyplot(tpr)
densityplot(tpr, main="densityplot( profile(lmer(..)) )")
splom(tpr)
\donttest{% for time constraint
## Batch and residual variance only
system.time(tpr2 <- profile(fm01ML, which=1:2, optimizer="Nelder_Mead"))
## GLMM example (running time ~11 seconds on a modern machine)
gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
            data = cbpp, family = binomial)
system.time(pr4 <- profile(gm1))
xyplot(pr4,layout=c(5,1),as.table=TRUE)
splom(pr4)
}%donttest
}
\keyword{methods}

