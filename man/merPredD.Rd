\name{merPredD}
\alias{merPredD}
\title{Generator object for the \code{\linkS4class{merPredD}} class}
\arguments{
  \item{X}{dense model matrix for the fixed-effects
  parameters, to be stored in the \code{X} field.}

  \item{Zt}{transpose of the sparse model matrix for the
  random effects.  It is stored in the \code{Zt} field.}

  \item{Lambdat}{transpose of the sparse lower triangular
  relative variance factor (stored in the \code{Lambdat}
  field).}

  \item{Lind}{integer vector of the same length as the
  \code{"x"} slot in the \code{Lambdat} field.  Its
  elements should be in the range 1 to the length of the
  \code{theta} field.}

  \item{theta}{numeric vector of variance component
  parameters (stored in the \code{theta} field).}
}
\description{
  The generator object for the
  \code{\linkS4class{merPredD}} reference class. Such an
  object is primarily used through its \code{new} method.
}
\section{Methods}{
  \describe{ \item{new(X, Zt, Lambdat, Lind,
  theta):}{Create a new \code{\linkS4class{merPredD}}
  object} }
}
\seealso{
  \code{\linkS4class{merPredD}}
}
\keyword{classes}

