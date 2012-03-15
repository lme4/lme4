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
  Substitute the '+' function for the '|' function in a
  mixed-model formula.  This provides a formula suitable
  for the current model.frame function.
}
\section{Note}{
  This function is called recursively on individual terms
  in the model, which is why the argument is called
  \code{term} and not a name like \code{form}, indicating a
  formula.
}
\examples{
subbars(Reaction ~ Days + (Days|Subject)) ## => Reaction ~ Days + (Days + Subject)
}
\seealso{
  \code{\link{formula}}, \code{\link{model.frame}},
  \code{\link{model.matrix}}.

  Other utilities: \code{\link{findbars}},
  \code{\link{mkRespMod}}, \code{\link{mkReTrms}},
  \code{\link{nlformula}}, \code{\link{nobars}}
}
\keyword{models}
\keyword{utilities}

