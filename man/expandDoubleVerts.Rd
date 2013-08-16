\name{expandDoubleVerts}
\alias{expandDoubleVerts}
\title{Expand terms with \code{'||'} notation into separate \code{'|'} terms}
\usage{
  expandDoubleVerts(term)
}
\arguments{
  \item{term}{a mixed-model formula}
}
\value{
 the modified term
}
\description{
 From the right hand side of a formula for a mixed-effects model,
 expand terms with the double vertical bar operator
 into separate, independent random effect terms.
}
\examples{
  f <- y ~ x + (x || g)
  # the right-hand side of f is,
  f[[3]]
  # the expanded right-hand side,
  expandDoubleVerts(f[[3]])
}
\seealso{
  \code{\link{formula}}, \code{\link{model.frame}},
  \code{\link{model.matrix}}.

  Other utilities: \code{\link{mkRespMod}},
  \code{\link{mkReTrms}}, \code{\link{nlformula}},
  \code{\link{nobars}}, \code{\link{subbars}}
}
\keyword{models}
\keyword{utilities}
