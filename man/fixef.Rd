\docType{methods}
\name{fixef}
\alias{fixed.effects}
\alias{fixef}
\alias{fixef.merMod}
\title{Extract fixed-effects estimates}
\usage{
  \method{fixef}{merMod} (object, ...)
}
\arguments{
  \item{object}{any fitted model object from which fixed
  effects estimates can be extracted.}

  \item{\dots}{optional additional arguments. Currently
  none are used in any methods.}
}
\value{
  a named, numeric vector of fixed-effects estimates.
}
\description{
  Extract the fixed-effects estimates
}
\details{
  Extract the estimates of the fixed-effects parameters
  from a fitted model.
}
\examples{
fixef(lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))
}
\keyword{models}

