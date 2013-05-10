\name{golden}
\alias{golden}
\title{Generator object for the golden search optimizer class.}
\usage{
  golden(...)
}
\arguments{
  \item{\dots}{additional, optional arguments.  None are
  used at present.}
}
\description{
  The generator objects for the \code{\linkS4class{golden}}
  class of a scalar optimizer for a parameter within an
  interval.  The optimizer uses reverse communications.
}
\note{
  Arguments to the \code{new} methods must be named
  arguments. \code{lower} and \code{upper} are the bounds
  for the scalar parameter; they must be finite.
}
\section{Methods}{
  \describe{ \item{\code{new(lower=lower,
  upper=upper)}}{Create a new \code{\linkS4class{golden}}
  object.} }
}
\seealso{
  \code{\linkS4class{golden}}
}
\keyword{classes}

