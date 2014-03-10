\name{ngrps}
\alias{ngrps}
\title{Number of levels}
\description{
  Number of levels of one or more grouping factor.
}
\usage{
ngrps(object, ...)
}
\arguments{
  \item{object}{An \code{R} object.}
  \item{...}{Currently ignored.}
}
\value{
  Number(s) of levels
}
\details{
  Currently there are methods for \code{merMod} and \code{factor}
  objects.
}
\examples{
ngrps(factor(seq(1,10,2)))
ngrps(lmer(Reaction ~ 1|Subject, sleepstudy))
}
