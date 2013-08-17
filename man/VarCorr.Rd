\name{VarCorr}
\title{Extract Variance and Correlation Components}
\alias{VarCorr}
\alias{VarCorr.merMod}
\usage{
\method{VarCorr}{merMod} (x, sigma, rdig)
}
\arguments{
  \item{x}{a fitted model object, usually an object inheriting from
    class \code{\linkS4class{merMod}}.}
  \item{sigma}{an optional numeric value used as a multiplier for the
    standard deviations.  Default is \code{1}.}
  \item{rdig}{an optional integer value specifying the number of digits
    used to represent correlation estimates.  Default is \code{3}.}
}
\description{
  This function calculates the estimated variances,
  standard deviations, and correlations between the
  random-effects terms in a mixed-effects model, of class
  \code{\linkS4class{merMod}} (linear, generalized or
  nonlinear).  The within-group error variance and standard
  deviation are also calculated.
}
\value{
  a list of matrices, one for each random effects grouping
  term.  For each grouping term, the standard deviations and
  correlation matrices for each grouping term are stored as
  attributes \code{"stddev"} and \code{"correlation"},
  respectively, of the variance-covariance matrix, and the
  residual standard deviation is stored as attribute
  \code{"sc"} (for \code{glmer} fits, this attribute stores
  the scale parameter of the model).
}
\author{
  This is modeled after \code{\link[nlme]{VarCorr}} from
  package \pkg{nlme}, by Jose Pinheiro and Douglas Bates.
}
\seealso{
  \code{\link{lmer}}, \code{\link{nlmer}}
}
\examples{
data(Orthodont, package="nlme")
fm1 <- lmer(distance ~ age + (age|Subject), data = Orthodont)
VarCorr(fm1)
}
\keyword{models}

