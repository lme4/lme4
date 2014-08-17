\name{splitForm}
\alias{splitForm}
\title{Split up a flexible mixed model formula into an unstructured and
  structured part}
\usage{
  splitForm(formula, specials = c("d", "cs", "ar1d"))
}
\arguments{
  \item{formula}{a mixed-model formula}
  \item{specials}{special functions for generating a structured
    random-effects term}
}
\description{
  Split up a flexible mixed model formula into an unstructured and
  structured part.
}
\value{
  A list with elements \code{lmerformula} containing the unstructured
  part of \code{formula} and \code{reGenerators} containing the
  structured part.
}
\seealso{
  \code{\link{formula}}, \code{\link{model.frame}},
  \code{\link{model.matrix}}.

  Other utilities: \code{\link{mkRespMod}},
  \code{\link{mkReTrms}}, \code{\link{nlformula}},
  \code{\link{nobars}}, \code{\link{subbars}}, \code{\link{findbars}}
}
\examples{
splitForm(f1 <- Reaction ~ Days + d(~(Days | Subject)), "d")
}
\keyword{models}
\keyword{utilities}
