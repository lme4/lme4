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
  \item{ord}{scalar integer between 1 and 100 - the order,
    or number of nodes and weights, in the rule.  When the
    function being multiplied by the standard normal density
    is a polynomial of order \eqn{2k-1} the rule of order \eqn{k}
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

  Originally based on package \CRANpkg{SparseGrid}'s hidden \code{GQN()},
  then on \CRANpkg{fastGHQuad}'s \code{gaussHermiteData(.)}.% ==> in ../R/GHrule.R
}
\references{
  Qing Liu and Donald A. Pierce (1994).
  A Note on Gauss-Hermite Quadrature. \emph{Biometrika} \bold{81}(3), 624--629.
  \doi{10.2307/2337136}
}
\seealso{a different interface is available via \code{\link{GQdk}()}.
}
\examples{
(r5  <- GHrule( 5, asMatrix=FALSE))
(r12 <- GHrule(12, asMatrix=FALSE))

## second, fourth, sixth, eighth and tenth central moments of the
## standard Gaussian N(0,1) density:
ps <- seq(2, 10, by = 2)
cbind(p = ps, "E[X^p]" = with(r5,  sapply(ps, function(p) sum(w * z^p)))) # p=10 is wrong for 5-rule
p <- 1:15
GQ12 <- with(r12, sapply(p, function(p) sum(w * z^p)))
cbind(p = p, "E[X^p]" = zapsmall(GQ12))
## standard R numerical integration can do it too:
intL <- lapply(p, function(p) integrate(function(x) x^p * dnorm(x),
                                        -Inf, Inf, rel.tol=1e-11))
integR <- sapply(intL, `[[`, "value")
cbind(p, "E[X^p]" = integR)# no zapsmall() needed here
all.equal(GQ12, integR, tol=0)# => shows small difference
stopifnot(all.equal(GQ12, integR, tol = 1e-10))
(xactMom <- cumprod(seq(1,13, by=2)))
stopifnot(all.equal(xactMom, GQ12[2*(1:7)], tol=1e-14))
## mean relative errors :
mean(abs(GQ12  [2*(1:7)] / xactMom - 1)) # 3.17e-16
mean(abs(integR[2*(1:7)] / xactMom - 1)) # 9.52e-17 {even better}
}

