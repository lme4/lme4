\name{mkdevfun}
\alias{mkdevfun}
\title{
  Create a deviance evaluation function from a merMod object
}
\description{
  From an merMod object create an R function that takes a single
  argument, which is the new parameter value, and returns the
  deviance.
}
\usage{
mkdevfun(pp, resp, rho, nAGQ=1L)
}
\arguments{
  \item{pp}{a prediction module, typically of class
    \code{\linkS4class{merPredD}}.%% FIXME merPred --> *D and *S ??
  }
  \item{resp}{a response module, e.g., of class
    \code{\linkS4class{lmerResp}}.}
  \item{rho}{an environment containing auxillary information for the
    function evaluation}
  \item{nAGQ}{scalar integer - the number of adaptive Gauss-Hermite
    quadrature points.  A value of 0 indicates that both the
    fixed-effects parameters and the random effects are optimized by the
    iteratively reweighted least squares algorithm.}
}
\details{
  The function returned by \code{mkdevfun} evaluates the deviance of the
  model represented by the predictor module, \code{pp}, and the response
  module, \code{resp}.

  For \code{\link{lmer}} model objects the argument of the resulting
  function is the variance component parameter, \code{theta}, with lower
  bound.  For \code{glmer} or \code{nlmer} model objects with
  \code{nAGQ = 0} the argument is also \code{theta}.  However, when {nAGQ > 0} the
  argument is \code{c(theta, beta)}.
}
\value{
  A function of one numeric argument.
}
%\references{}
%\author{}
%\note{}
\seealso{\code{\link{lmer}}, \code{\link{glmer}} and \code{\link{nlmer}}
}
\examples{
(dd <- lmer(Yield ~ 1|Batch, Dyestuff, devFunOnly=TRUE))
dd(0.8)
bobyqa(1, dd, 0)
}
\keyword{models}
