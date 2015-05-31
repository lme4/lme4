\docType{methods}
\name{fixef}
\alias{fixed.effects}
\alias{fixef}
\alias{fixef.merMod}
\title{Extract fixed-effects estimates}
\usage{
  \method{fixef}{merMod} (object, add.dropped=FALSE, ...)
}
\arguments{
  \item{object}{any fitted model object from which fixed
  effects estimates can be extracted.}
\item{add.dropped}{for models with rank-deficient design
  matrix, reconstitute the full-length parameter vector by
  adding \code{NA} values in appropriate locations?}
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
fm2 <- lmer(Reaction ~ Days + Days2 + (1|Subject),
            data=transform(sleepstudy,Days2=Days))
fixef(fm2,add.dropped=TRUE)
## first two parameters are the same ...
stopifnot(all.equal(fixef(fm2,add.dropped=TRUE)[1:2],
                    fixef(fm2)))
}
\keyword{models}

