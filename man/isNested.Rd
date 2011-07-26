\name{isNested}
\alias{isNested}
\title{Is f1 nested within f2?}
\usage{isNested(f1, f2)}
\description{Is f1 nested within f2?}
\arguments{
  \item{f1}{\code{\link{factor}} 1}
  \item{f2}{factor 2}
}
\details{Does every level of \code{f1} occur in conjunction with exactly one
  level of \code{f2}? The function is based on converting a triplet sparse
  matrix to a compressed column-oriented form in which the nesting can
  be quickly evaluated.}
\value{\code{TRUE} if factor 1 is nested within factor 2, else \code{FALSE}}
\examples{
with(Pastes, isNested(sample, batch))  ## => TRUE
with(Pastes, isNested(cask, batch))    ## => FALSE
}
\keyword{manip}
