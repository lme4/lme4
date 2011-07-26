\name{mkdevfun}
\alias{mkdevfun}
\title{
  Create a deviance evaluation function from an merMod object
}
\description{
  From an merMod object create an R function that takes a single
  argument, which is the new parameter value, and returns the
  deviance.
}
\usage{
mkdevfun(mod, nAGQ = 1L, beta0 = numeric(length(mod@fe@coef)),
         u0 = numeric(length(mod@re@u)),
         verbose = 0L, compDev = TRUE)
}
\arguments{
  \item{mod}{
    an object that inherits from class \code{\linkS4class{merMod}}.
  }
  \item{nAGQ}{
    number of points per axis for adaptive Gauss-Hermite quadrature. The
    default, \code{1}, creates the Laplace approximation.  The special
    value, \code{0L}, creates the PIRLSBeta Laplace approximation in
    which the penalized iteratively reweighted least squares (PIRLS) algorithm
    is used to optimize all the coefficients.
  }
  \item{beta0}{
    Starting value to use for \code{beta} in the PIRLS or PIRLSBeta algorithm
  }
  \item{u0}{
    Starting value to use for \code{u} in the PIRLS or PIRLSBeta algorithm
  }
  \item{verbose}{
    integer value indicating level of verbosity in the output.  Values
    greater than 1 cause iteration output during the PIRLS algorithm,
    which can be very verbose.
  }
  \item{compDev}{
    should the compiled
    deviance evaluation be used?  A value of \code{FALSE} provides an
    evaluation function using functions from the Matrix package only.
    (Only available for \code{lmer} model objects).
}
}
\details{
  The function returned by \code{mkdevfun} evaluates the deviance but
  does not update any slots in \code{mod}, thus it can be used to fit
  variations on a fitted model without modifying the original object.
  
  For \code{lmer} model objects the argument of the resulting
  function is the variance component parameter, \code{theta}, with lower
  bound \code{mod@re@lower}.  For \code{glmer} or
  \code{nlmer} model objects with \code{nAGQ = 0} the
  argument is also \code{theta}.  However, when {nAGQ > 0} the argument
  is \code{c(theta, beta)}.
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
fm1 <- lmer(Yield ~ 1|Batch, Dyestuff, doFit = FALSE)
(dd0 <- mkdevfun(fm1))
dd0(0.8)
bobyqa(1, dd0, 0)
(dd1 <- mkdevfun(fm1, compDev = FALSE))
\dontrun{
dd1(0.8)
bobyqa(1, dd1, 0)
}
}
\keyword{models}
