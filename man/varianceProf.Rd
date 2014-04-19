\name{varianceProf}
\alias{varianceProf}
\title{Transform Profile to the variance scale}
\description{
  Transform a mixed-effects profile, i.e. conceptually,
  \code{x <- \link{profile}(\link{lmer}(....),  ..)}
  from standard deviation to the variance scale.
}
\usage{
varianceProf(x, ranef = TRUE)
}
\arguments{
  \item{x}{a mixed-effects model profile.}
  \item{ranef}{logical indicating if all sigmas should be transformed,
    see \code{\link{logProf}}.}
}
\value{
  a transformed mixed-effects model profile
}
\seealso{
  \code{\link[=profile.merMod]{profile}} and methods, notably,
  \code{\link{logProf}} for log transformation of the sigmas;
  also for nice examples including graphics.
}
