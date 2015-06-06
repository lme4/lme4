\name{merPredD}
\alias{merPredD}
\title{Generator object for the \code{\linkS4class{merPredD}} class}
\usage{
merPredD(...)
}
\arguments{
  \item{...}{List of arguments (see Note).}
}
\description{
  The generator object for the
  \code{\linkS4class{merPredD}} reference class. Such an
  object is primarily used through its \code{new} method.
}
\note{
  \code{merPredD(...)} is a short form of \code{new("merPredD", ...)}
  to create a new \code{\linkS4class{merPredD}} object and the
  \code{...} must be named arguments, \code{(X, Zt, Lambdat, Lind,
    theta,n)}:
  \describe{
    \item{\code{X}:}{dense model matrix for the fixed-effects parameters, to be
      stored in the \code{X} field.}
    \item{\code{Zt}:}{transpose of the sparse model matrix for the random
      effects.  It is stored in the \code{Zt} field.}
    \item{\code{Lambdat}:}{transpose of the sparse lower triangular relative
      variance factor (stored in the \code{Lambdat} field).}
    \item{\code{Lind}:}{integer vector of the same length as the \code{x} slot in the
      \code{Lambdat} field.  Its elements should be in the range 1 to
      the length of the \code{theta} field.}
    \item{\code{theta}:}{numeric vector of variance component parameters (stored
      in the \code{theta} field).}
    \item{\code{n}:}{sample size, usually \code{nrow(X)}.}
  }
}
\seealso{
  The class definition, \code{\linkS4class{merPredD}}, also for examples.
}
\keyword{classes}

