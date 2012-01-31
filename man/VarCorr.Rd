\name{VarCorr}
\alias{VarCorr}
\alias{VarCorr.merMod}
\title{Extract variance and correlation components}
\usage{
  \method{VarCorr}{merMod} (x, sigma, rdig)
}
\arguments{
  \item{x}{a fitted model object, usually an object
  inheriting from class \code{\linkS4class{merMod}}.}

  \item{sigma}{an optional numeric value used as a
  multiplier for the standard deviations.  Default is
  \code{1}.}

  \item{rdig}{an optional integer value specifying the
  number of digits used to represent correlation estimates.
  Default is \code{3}.}
}
\value{
  a matrix with the estimated variances, standard
  deviations, and correlations for the random effects.  The
  first two columns, named \code{Variance} and
  \code{StdDev}, give, respectively, the variance and the
  standard deviations.  If there are correlation components
  in the random effects model, the third column, named
  \code{Corr}, and the remaining unnamed columns give the
  estimated correlations among random effects within the
  same level of grouping.  The within-group error variance
  and standard deviation are included as the last row in
  the matrix.
}
\description{
  This function calculates the estimated variances,
  standard deviations, and correlations between the
  random-effects terms in a mixed-effects model, of class
  \code{\linkS4class{merMod}} (linear, generalized or
  nonlinear).  The within-group error variance and standard
  deviation are also calculated.
}
\examples{
data(Orthodont, package="nlme")
fm1 <- lmer(distance ~ age + (age|Subject), data = Orthodont)
VarCorr(fm1)
}
\author{
  This is modeled after \code{\link[nlme]{VarCorr}} from
  package \pkg{nlme}, by Jose Pinheiro and Douglas Bates.
}
\seealso{
  \code{\link{lmer}}, \code{\link{nlmer}}
}
\keyword{models}

