\name{GHrule}
\alias{GHrule}
\title{
  Nodes and weights for Gauss-Hermite quadrature
}
\description{
  Returns the nodes, weights and the logarithm of the standard Gaussian
  density for a Gauss-Hermite quadrature rule, with the standard
  Gaussian density as the kernel, of order \code{ord}.
}
\usage{
GHrule(ord)
}
\arguments{
  \item{ord}{
    positive integer not exceeding 25 - the order of the rule.
  }
}
\details{
  Univariate Gauss-Hermite quadrature evaluates the integral of a
  function that is multiplied by a \dQuote{kernel} where the kernel is a
  multiple of \eqn{e^{-z^2}}{exp(-z^2)} or
  \eqn{e^{-z^2/2}}{exp(-z^2/2)}.  For statisticians the natural
  candidate is the standard normal density,
  \eqn{\phi(z)=e^{-z^2/2}/\sqrt(2\pi)}{exp(-z^2/2)/sqrt(2*pi)}.
  A \eqn{k}{k}th-order Gauss-Hermite formula provides knots,
  \eqn{z_i,i=1,...,k}{z[i]}, and weights, \eqn{w_i,i=1,\dots,k}{w[i]}, such that
  \deqn{\int_{\mathbb{R}}t(z)\phi(z)\,dz\approx\sum_{i=1}^kw_it(z_i)}{int[-Inf]^Inf t(z) phi(z) dz ~ sum[i=1]^k w[i] t(z[i])}
}
\note{
  The choice of the value of \eqn{k} depends on the behavior of the function
  \eqn{t(z)}.  If \eqn{t(z)} is a polynomial of degree \eqn{k-1} then the
  Gauss-Hermite formula for order \eqn{k} or greater provides the exact
  value of the integral.  The fact that we want \eqn{t(z)} to behave
  like a low-order polynomial is often neglected in the formulation of
  a Gauss-Hermite approximation to an integral.
}
\value{
  a numeric matrix with \code{ord} rows and three columns named
  \dQuote{z}, \dQuote{w} and \dQuote{ldnorm}.
}
\references{
  Based on code in the \dQuote{SparseGrid} package.
}
%\author{}
%\note{}
%\seealso{}
\examples{
GHrule(5)
}
\keyword{math}
