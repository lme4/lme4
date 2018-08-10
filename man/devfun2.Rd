\name{devfun2}
\alias{devfun2}
\title{
  Deviance function defined in terms of standard deviations/correlations
}
\description{
The deviance is profiled with respect to the fixed-effects
parameters but not with respect to sigma; that is, the
function takes parameters for the variance-covariance parameters
and for the residual standard deviation. The random-effects
variance-covariance parameters
are on the standard deviation/correlation scale, not the theta (Cholesky factor)
scale.
}
\usage{
devfun2(fm, useSc,
   transfuns = list(from.chol = Cv_to_Sv,
                    to.chol = Sv_to_Cv, to.sd = identity), ...)
}
\arguments{
  \item{fm}{a fitted model of class \sQuote{merMod}}
  \item{useSc}{
    (logical) whether a scale parameter is used
  }
  \item{transfuns}{
    functions for converting parameters to and from the
    Cholesky-factor scale
  }
  \item{\dots}{
    arguments passed to the internal \code{profnames} function
    (\code{signames=TRUE} to use old-style \code{.sigxx names},
    FALSE uses (sd_cor|xx);
    also \code{prefix=c("sd","cor")})
  }
}
\value{
  Returns a function that takes a vector of standard deviations and
  correlations and returns the deviance (or REML criterion). The
  function has additional attributes
  \describe{
    \item{optimum}{a named vector giving the parameter values
      at the optimum}
    \item{basedev}{the deviance at the optimum (\emph{not} the
      REML criterion, even if the original model was fitted using
      \code{REML=TRUE})}
    \item{thopt}{the optimal variance-covariance parameters on the
      \dQuote{theta} (Cholesky factor) scale}
    \item{stderr}{standard errors of fixed effect parameters}
  }
}
\examples{
m1 <- lmer(Reaction~Days+(Days|Subject),sleepstudy)
dd <- devfun2(m1,useSc=TRUE)
pp <- attr(dd,"optimum")
## extract variance-covariance and residual std dev parameters
sigpars <- pp[grepl("^\\\\.sig",names(pp))]
all.equal(unname(dd(sigpars)),deviance(refitML(m1)))
}
\keyword{utilities}
