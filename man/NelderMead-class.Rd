\docType{class}
\name{NelderMead-class}
\alias{NelderMead}
\alias{NelderMead-class}
\title{Class \code{"NelderMead"} of Nelder-Mead optimizers and its Generator}
\description{
  Class \code{"NelderMead"} is a reference class for a Nelder-Mead
  simplex optimizer allowing box constraints on the parameters and using
  reverse communication.

  The \code{NelderMead()} function conveniently generates such objects.
}
\usage{
NelderMead(...)
}
\arguments{
  \item{\dots}{Argument list (see Note below).}
}
\note{
  This is the default optimizer for the second stage of
  \code{\link{glmer}} and \code{\link{nlmer}} fits.  We
  found that it was more reliable and often faster than
  more sophisticated optimizers.

  Arguments to \code{NelderMead()} and the \code{\link{new}} method must
  be named arguments:
  \describe{
    \item{lower}{numeric vector of lower bounds - elements may be \code{-Inf}.}
    \item{upper}{numeric vector of upper bounds - elements may be \code{Inf}.}
    \item{xst}{numeric vector of initial step sizes to establish the
      simplex - all elements must be non-zero.}
    \item{x0}{numeric vector of starting values for the parameters.}
    \item{xt}{numeric vector of tolerances on the parameters.}
  }
}
\section{Methods}{
  \describe{\code{NelderMead$new(lower, upper, xst, x0, xt)}}{Create a
    new \code{\linkS4class{NelderMead}} object}
}
\section{Extends}{
  All reference classes extend and inherit methods from
  \code{"\linkS4class{envRefClass}"}.
}
\references{
  Based on code in the NLopt collection.
}
\seealso{
  \code{\link{Nelder_Mead}}, the typical \dQuote{constructor}.
  Further, \code{\link{glmer}}, \code{\link{nlmer}}
}
\examples{
showClass("NelderMead")
}
\keyword{classes}

