\name{NelderMead}
\alias{NelderMead}
\title{Generator object for the Nelder-Mead optimizer class.}
\arguments{
  \item{lower}{numeric vector of lower bounds - elements
  may be \code{-Inf}.}

  \item{upper}{numeric vector of upper bounds - elements
  may be \code{Inf}.}

  \item{xst}{numeric vector of initial step sizes to
  establish the simplex - all elements must be non-zero.}

  \item{x0}{numeric vector of starting values for the
  parameters.}

  \item{xt}{numeric vector of tolerances on the
  parameters.}

  \item{\dots}{additional, optional arguments.  None are
  used at present.}
}
\description{
  The generator objects for the
  \code{\linkS4class{NelderMead}} class of optimizers
  subject to box constraints and using reverse
  communications.
}
\note{
  Arguments to the \code{new} methods must be named
  arguments.
}
\section{Methods}{
  \describe{\code{new(lower, upper, xst, x0, xt)}}{Create a
  new \code{\linkS4class{NelderMead}} object}
}
\seealso{
  \code{\linkS4class{NelderMead}}
}
\keyword{classes}

