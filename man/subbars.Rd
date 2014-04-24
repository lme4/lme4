\name{subbars}
\alias{subbars}
\title{"Sub[stitute] Bars"}
\usage{
  subbars(term)
}
\arguments{
  \item{term}{a mixed-model formula}
}
\value{
  the formula with all | operators replaced by +
}
\description{
  Substitute the '+' function for the '|' function in a mixed-model
  formula, recursively (hence the argument name \code{term}).  This
  provides a formula suitable for the current \code{\link{model.frame}}
  function.
}
\seealso{
  \code{\link{formula}}, \code{\link{model.frame}},
  \code{\link{model.matrix}}.

  Other utilities: \code{\link{findbars}}, \code{\link{nobars}},
  \code{\link{mkRespMod}}, \code{\link{mkReTrms}},
  \code{\link{nlformula}}.
}
\examples{
subbars(Reaction ~ Days + (Days|Subject)) ## => Reaction ~ Days + (Days + Subject)
}
\keyword{models}
\keyword{utilities}

