\docType{class}
\name{NelderMead-class}
\alias{NelderMead-class}
\title{Class \code{"NelderMead"}}
\description{
  A reference class for a Nelder-Mead simplex optimizer
  allowing box constraints on the parameters and using
  reverse communication.
}
\note{
  This is the default optimizer for the second stage of
  \code{\link{glmer}} and \code{\link{nlmer}} fits.  We
  found that it was more reliable and often faster than
  more sophisticated optimizers.
}
\section{Extends}{
  All reference classes extend and inherit methods from
  \code{"\linkS4class{envRefClass}"}.
}
\examples{
showClass("NelderMead")
}
\references{
  Based on code in the NLopt collection.
}
\seealso{
  \code{\link{glmer}}, \code{\link{nlmer}}
}
\keyword{classes}

