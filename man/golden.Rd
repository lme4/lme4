\name{golden}
\alias{golden}
\title{Generator object for the Nelder-Mead optimizer class.}
\section{Methods}{
  \describe{
    \item{\code{new(lower, upper, \dots)}:}{Create a new
      \code{\linkS4class{golden}} object.}
  }
}
\arguments{
  \item{lower}{lower bound for the scalar parameter - must be finite.}
  \item{upper}{upper bound for the scalar parameter - must be finite.}
  \item{\dots}{additional, optional arguments.  None are used at present.}
}
\note{Arguments to the \code{new} methods must be named arguments.}
\description{The generator objects for the
  \code{\linkS4class{golden}} class of a scalar optimizer for a
  parameter within an interval.  The optimizer uses reverse communications.}
\seealso{\code{\linkS4class{golden}}}
\keyword{classes}
