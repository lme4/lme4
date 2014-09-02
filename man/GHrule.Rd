\name{GHrule}
\alias{GHrule}
\title{Univariate Gauss-Hermite quadrature rule}
\description{
  Create a univariate Gauss-Hermite quadrature rule.
}
\usage{
  GHrule(ord, asMatrix = TRUE)
}
\arguments{
  \item{ord}{scalar integer between 1 and 25 - the order,
    or number of nodes and weights, in the rule.  When the
    function being multiplied by the standard normal density
    is a polynomial of order 2k-1 the rule of order k
    integrates the product exactly.}

  \item{asMatrix}{logical scalar - should the result be
    returned as a matrix.  If \code{FALSE} a data frame is
    returned.  Defaults to \code{TRUE}.}
}
\value{
  a matrix (or data frame, is \code{asMatrix} is false) with \code{ord}
  rows and three columns which are \code{z} the node positions, \code{w}
  the weights and \code{ldnorm}, the logarithm of the normal density
  evaluated at the nodes.
}
\details{
  This version of Gauss-Hermite quadrature provides the node positions
  and weights for a scalar integral of a function multiplied by the
  standard normal density.

  Originally based on package \pkg{SparseGrid}'s hidden \code{GQN()}.
}
\seealso{a different interface is available via \code{\link{GQdk}()}.
}
\examples{
(r5 <- GHrule(5, asMatrix=FALSE))
## second, fourth, sixth, eighth and tenth central moments of the
## standard Gaussian density
with(r5, sapply(seq(2, 10, 2), function(p) sum(w * z^p)))
}

