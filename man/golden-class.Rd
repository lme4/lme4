\name{golden-class}
\Rdversion{1.1}
\docType{class}
\alias{golden-class}
\title{Class \code{"golden"}}
\description{
  A reference class for a golden search scalar optimizer using reverse
  communication.
}
\section{Extends}{
  All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.
}
%\references{}
%\author{}
%\note{}
%\seealso{}
\examples{
showClass("golden")
}
\keyword{classes}
\section{Fields}{
  \describe{
    \item{\code{Ptr}:}{\code{externalptr} to the instance of the C++ class}
    \item{\code{lowerbd}:}{numeric scalar - the lower bound on the
      interval within which to minimize the objective.}
    \item{\code{upperbd}:}{numeric scalar - the upper bound on the
      interval within which to minimize the objective.}
  }
}
\section{Methods}{
  \describe{
    \item{\code{newf(value)}:}{provide the function value at the current
      \code{xeval}.}
    \item{\code{value()}:}{return the lowest function value recorded to
      this point.}
    \item{\code{ptr()}:}{return the external pointer, \code{Ptr},
      regenerating the C++ object if necessary.}
    \item{\code{xeval()}:}{next value of the parameter at which to
      evaluate the objective function.}
    \item{\code{xpos()}:}{position of the minimum objective value seen
      to this point.}
  }
}
