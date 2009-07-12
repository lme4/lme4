\name{fixef}
\docType{genericFunction}
\alias{fixef}
\alias{fixed.effects}
\alias{fixef,ANY-method}
\alias{fixef,mer-method}
\title{Extract the fixed-effects estimates}
\description{
  Extract the estimates of the fixed-effects parameters from a fitted model.
}
\usage{
fixef(object, \dots)
fixed.effects(object, \dots)  # deprecated
}
\arguments{
 \item{object}{any fitted model object from which fixed effects
   estimates can be extracted.}
 \item{\dots}{optional additional arguments. Currently none are used in
   any methods.}
}
\value{
  a named, numeric vector of fixed-effects estimates.
}
\examples{
fixef(lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))
}
\keyword{models}
