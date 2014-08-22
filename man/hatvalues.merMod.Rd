\name{hatvalues.merMod}
\alias{hatvalues.merMod}
\title{Diagonal elements of the hat matrix}
\usage{
  \method{hatvalues}{merMod} (model, fullHatMatrix = FALSE, ...)
}
\arguments{
  \item{model}{An object of class \code{\link{merMod}}.}
  \item{fullHatMatrix}{Return full hat matrix (not just diagonal values)?}
  \item{...}{Not currently used}
}
\value{
  The diagonal elements of the hat matrix.
}
\description{
  Returns the values on the diagonal of the hat matrix, which is the
  matrix that transforms the response vector (minus any offset) into the
  fitted values (minus any offset).  Note that this method should only
  be used for linear mixed models.  It is not clear if the hat matrix
  concept even makes sense for generalized linear mixed models.
}
\examples{
m <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
hatvalues(m)
}
