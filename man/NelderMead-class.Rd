\name{NelderMead-class}
\Rdversion{1.1}
\docType{class}
\alias{NelderMead-class}
\title{Class \code{"NelderMead"}}
\description{
  A reference class for a Nelder-Mead simplex optimizer allowing box
  constraints on the parameters and using reverse communication.
}
\section{Extends}{
  All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.
}
\references{
  Based on code in the NLopt collection.
}
%\author{}
\note{This is the default optimizer for the second stage of
  \code{\link{glmer}} and \code{\link{nlmer}} fits.  We found that it
  was more reliable and often faster than more sophisticated optimizers.}
\seealso{\code{\link{glmer}}, \code{\link{nlmer}}}
\examples{
showClass("NelderMead")
}
\keyword{classes}
\section{Fields}{
  \describe{
    \item{\code{Ptr}:}{\code{externalptr} to the instance of the C++ class}
    \item{\code{lowerbd}:}{\code{numeric} vector of lower bounds - some
      or all of the elements may be \code{-Inf} for parameters with no
      lower bound.}
    \item{\code{upperbd}:}{\code{numeric} vector of upper bounds - some
      or all of the elements may be \code{Inf} for parameters with no
      upper bound.}
    \item{\code{xstep}:}{\code{numeric} vector of initial steps to
      generate the simplex.  All the elements must be non-zero.}
    \item{\code{xtol}:}{\code{numeric} vector of tolerances on the parameters}
  }
}
\section{Methods}{
  \describe{
    \item{\code{setFtolAbs(fta)}:}{set the absolute tolerance on the
      function value.}
    \item{\code{newf(value)}:}{provide the function value at the current
      \code{xeval} vector.}
    \item{\code{setMaxeval(mxev)}:}{set the maximum number of function evaluations.}
    \item{\code{value()}:}{return the lowest function value previously recorded}
    \item{\code{setMinfMax(minf)}:}{set the maximum value to be declared
      as the minimum function value.  Defaults to \code{-Inf}.}
    \item{\code{setFtolRel(ftr)}:}{set the relative tolerance on the
      function value.}
    \item{\code{setForceStop(stp)}:}{Force the iterations to stop.  This
      method is unlikely to be used in an algorithm using reverse communication.}
    \item{\code{ptr()}:}{return the external pointer, \code{Ptr},
      regenerating the C++ object if necessary.}
    \item{\code{xeval()}:}{next value of the parameters at which to
      evaluate the objective function}
    \item{\code{setIprint(iprint)}:}{set the parameter controlling
      printing of intermediate results.  If \code{iprint > 0} then the
      function and parameter values are printed every \code{iprint} evaluations.}
    \item{\code{xpos()}:}{position of the minimum objective value seen
      to this point.}
  }
}
