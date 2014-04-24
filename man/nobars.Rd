\name{nobars}
\alias{nobars}
\title{Omit terms separated by vertical bars in a formula}
\usage{
  nobars(term)
}
\arguments{
  \item{term}{the right-hand side of a mixed-model formula}
}
\value{
  the fixed-effects part of the formula
}
\description{
  Remove the random-effects terms from a mixed-effects
  formula, thereby producing the fixed-effects formula.
}
\note{
  This function is called recursively on individual terms
  in the model, which is why the argument is called
  \code{term} and not a name like \code{form}, indicating a
  formula.
}
\examples{
nobars(Reaction ~ Days + (Days|Subject)) ## => Reaction ~ Days
}
\seealso{
  \code{\link{formula}}, \code{\link{model.frame}},
  \code{\link{model.matrix}}.

  Other utilities: \code{\link{findbars}},
  \code{\link{mkRespMod}}, \code{\link{mkReTrms}},
  \code{\link{nlformula}}, \code{\link{subbars}}
}
\keyword{models}
\keyword{utilities}

