\name{profile-methods}
\docType{methods}
\alias{profile-methods}
\alias{profile,ANY-method}
\alias{profile,merMod-method}
\title{Methods for profile() of [ng]lmer fitted models}
\description{
  Methods for function \code{\link{profile}} (package \pkg{stats}),
  here for profiling (fitted) mixed effect models.
}
\usage{
\S4method{profile}{merMod}(fitted, alphamax = 0.01, maxpts = 100,
          delta = cutoff/8, tr = 0, \dots)
}
\arguments{
  \item{fitted}{a fitted model, e.g., the result of \code{\link{lmer}(..)}.}
  \item{alphamax}{used when \code{delta} is unspecified, as probability
    ... to compute \code{delta} ...}
  \item{maxpts}{ ... % FIXME
}
  \item{delta}{ ... % FIXME
}
  \item{tr}{ ... % FIXME
}
  \item{\dots}{potential further arguments for \code{profile} methods.}
}

\section{Methods}{
  FIXME: These (signatures) will change soon --- document \bold{after} change!

  \describe{
    \item{\code{signature(fitted = "merMod")}}{ ...
    }
  }
}
\seealso{
  For (more expensive) alternative confidence intervals:
  \code{\link{bootMer}}.
}
\examples{
\dontrun{
%% Do keep at least *one* such example! -- this is (also) a regression test
fm01ML <- lmer(Yield ~ 1|Batch, Dyestuff, REML = FALSE)

## 0.8s (on a 5600 MIPS 64bit fast(year 2009) desktop "AMD Phenom(tm) II X4 925"):
system.time( tpr <- profile(fm01ML) )

(confint(tpr) -> CIpr)
print(xyplot(tpr))
stopifnot(dim(CIpr) == c(3,2),
          all.equal(unname(CIpr[2,]),c(3.64362, 4.21446), tol=1e-6))
}
}
\keyword{methods}
