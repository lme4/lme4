\name{mkdevfun}
\alias{mkdevfun}
\title{Create Deviance Evaluation Function from Predictor and Response Module}
\description{
  From a \code{\linkS4class{merMod}} object create an \R function that
  takes a   single argument, which is the new parameter value, and
  returns the deviance.
}
\usage{
mkdevfun(rho, nAGQ = 1L, maxit = 100, verbose = 0, control = list())
}
\arguments{
  \item{rho}{an environment containing \code{pp}, a
    prediction module, typically of class
    \code{\linkS4class{merPredD}} and \code{resp}, a response
    module, e.g., of class \code{\linkS4class{lmerResp}}.}

  \item{nAGQ}{scalar integer - the number of adaptive
    Gauss-Hermite quadrature points.  A value of 0 indicates
    that both the fixed-effects parameters and the random
    effects are optimized by the iteratively reweighted least
    squares algorithm.}

  \item{maxit}{scalar integer, currently only for GLMMs: the maximal
    number of Pwrss update iterations.}
  \item{verbose}{scalar logical: print verbose output?}

  \item{control}{list of control parameters, a subset of
    those specified by \code{\link{lmerControl}}
    (\code{tolPwrss} and \code{compDev} for GLMMs,
    \code{tolPwrss} for NLMMs)}
}
\value{
  A \code{\link{function}} of one numeric argument.
}
\details{
  The function returned by \code{mkdevfun} evaluates the
  deviance of the model represented by the predictor
  module, \code{pp}, and the response module, \code{resp}.

  For \code{\link{lmer}} model objects the argument of the
  resulting function is the variance component parameter,
  \code{theta}, with lower bound.  For \code{\link{glmer}} or
  \code{\link{nlmer}} model objects with \code{nAGQ = 0} the
  argument is also \code{theta}.  However, when nAGQ > 0,
  the argument is \code{c(theta, beta)}.
}
\examples{
(dd <- lmer(Yield ~ 1|Batch, Dyestuff, devFunOnly=TRUE))
dd(0.8)
minqa::bobyqa(1, dd, 0)
}
\seealso{
  \code{\link{lmer}}, \code{\link{glmer}} and
  \code{\link{nlmer}}
}
\keyword{models}
