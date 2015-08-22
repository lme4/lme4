\name{golden-class}
\docType{class}
\title{Class \code{"golden"} and Generator for Golden Search Optimizer Class}
\alias{golden-class}
\alias{golden}
\description{
  \code{"golden"} is a reference class for a golden search scalar optimizer,
  for a parameter within an interval.

  \code{golden()} is the generator for the \code{"golden"}
  class.  The optimizer uses reverse communications.
}
\usage{
golden(...)
}
\arguments{
  \item{\dots}{(partly optional) arguments passed to
    \code{\link{new}()} must be named arguments.  \code{lower} and
    \code{upper} are the bounds for the scalar parameter; they must be
    finite.}
}
\section{Extends}{
  All reference classes extend and inherit methods from
  \code{"\linkS4class{envRefClass}"}.
}
\examples{
showClass("golden")

golden(lower= -100, upper= 1e100)
}
\keyword{classes}

