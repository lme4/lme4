\name{dummy}
\alias{dummy}
\title{Dummy variables (experimental)}
\usage{
  dummy(f, levelsToKeep)
}
\arguments{
  \item{f}{An object coercible to \code{\link{factor}}.}
  \item{levelsToKeep}{An optional character vector giving the subset of
  \code{levels(f)} to be converted to dummy variables.}
}
\value{
  A \code{\link{model.matrix}} with dummy variables as columns.
}
\description{
  Largely a wrapper for \code{\link{model.matrix}} that
  accepts a factor, \code{f}, and returns a dummy
  matrix with \code{nlevels(f)-1} columns (the first
  column is dropped by default).  Useful whenever one wishes to avoid
  the behaviour of \code{model.matrix} of always returning
  an \code{nlevels(f)}-column matrix, either by the
  addition of an intercept column, or by keeping one
  column for all levels.
}
\examples{
data(Orthodont,package="nlme")
lmer(distance ~ age + (age|Subject) +
     (0+dummy(Sex, "Female")|Subject), data = Orthodont)
}
