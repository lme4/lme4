\name{findbars}
\alias{findbars}
\title{Determine random-effects expressions from a formula}
\usage{
  findbars(term)
}
\arguments{
  \item{term}{a mixed-model formula}
}
\description{
  From the right hand side of a formula for a mixed-effects
  model, determine the pairs of expressions that are
  separated by the vertical bar operator.  Also expand the
  slash operator in grouping factor expressions and expand
  terms with the double vertical bar operator into separate,
  independent random effect terms.
}
\value{
  pairs of expressions that were separated by vertical bars
}
\note{
  This function is called recursively on individual terms
  in the model, which is why the argument is called
  \code{term} and not a name like \code{form}, indicating a
  formula.
}
\seealso{
  \code{\link{formula}}, \code{\link{model.frame}},
  \code{\link{model.matrix}}.

  Other utilities: \code{\link{mkRespMod}},
  \code{\link{mkReTrms}}, \code{\link{nlformula}},
  \code{\link{nobars}}, \code{\link{subbars}}
}
\examples{
findbars(f1 <- Reaction ~ Days + (Days | Subject))
## => list( Days | Subject )
## These two are equivalent:% tests in ../inst/tests/test-doubleVertNotation.R
findbars(y ~ Days + (1 | Subject) + (0 + Days | Subject))
findbars(y ~ Days + (Days || Subject))
## => list of length 2:  list ( 1 | Subject ,  0 + Days | Subject)
findbars(~ 1 + (1 | batch / cask))
## => list of length 2:  list ( 1 | cask:batch ,  1 | batch)
\dontshow{
stopifnot(identical(findbars(f1),
                    list(quote(Days | Subject))))
}
}
\keyword{models}
\keyword{utilities}

