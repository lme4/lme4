\name{rePos}
\alias{rePos}
\title{Generator object for the rePos (random-effects positions) class}
\usage{
  rePos(...)
}
\arguments{
  \item{\dots}{Argument list (see Note).}
}
\description{
  The generator object for the \code{\linkS4class{rePos}}
  class used to determine the positions and orders of
  random effects associated with particular random-effects
  terms in the model.
}
\note{
  Arguments to the \code{new} methods must be named
  arguments. \code{mer}, an object of class
  \code{"\linkS4class{merMod}"}, is the only
  required/expected argument.
}
\section{Methods}{
  \describe{ \item{\code{new(mer=mer)}}{Create a new
  \code{\linkS4class{rePos}} object.} }
}
\seealso{
  \code{\linkS4class{rePos}}
}
\keyword{classes}

