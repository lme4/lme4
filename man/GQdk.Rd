\name{GQdk}
\title{Sparse Gaussian / Gauss-Hermite Quadrature grid}
\alias{GQdk}
\alias{GQN}
\description{
  Generate the sparse multidimensional Gaussian quadrature grids.

  Currently unused.  See \code{\link{GHrule}()} for the version
  currently in use in package \pkg{lme4}.
}
\usage{
  GQdk(d = 1L, k = 1L)
  GQN
}
\arguments{
  \item{d}{integer scalar - the dimension of the function
    to be integrated with respect to the standard
    \code{d}-dimensional Gaussian density.}
  \item{k}{integer scalar - the order of the grid.  A grid
    of order \code{k} provides an exact result for a
    polynomial of total order of \code{2k - 1} or less.}
}
\value{
  \code{GQdk()} returns a matrix with \code{d + 1} columns.  The first
  column is the weights and the remaining \code{d} columns are the
  node coordinates.

  \code{GQN} is a \code{\link{list}} of lists, containing the
  non-redundant quadrature nodes and weights for integration of a scalar
  function of a \code{d}-dimensional argument with respect to the density
  function of the \code{d}-dimensional Gaussian density function.
  \cr
  The outer list is indexed by the dimension, \code{d}, in the
  range of 1 to 20.  The inner list is indexed by \code{k},
  the order of the quadrature.
}
\note{
  \code{GQN} contains only the non-redundant nodes.  To regenerate
  the whole array of nodes, all possible permutations of
  axes and all possible combinations of \eqn{\pm 1}{+/- 1}
  must be applied to the axes.  This entire array of nodes is exactly
  what \code{\link{GQdk}()} reproduces.

  The number of nodes gets very large very quickly with
  increasing \code{d} and \code{k}.  See the charts at
  \url{http://www.sparse-grids.de}.
}
\examples{
GQdk(2,5) # 53 x 3

GQN[[3]][[5]] # a 14 x 4 matrix
}

