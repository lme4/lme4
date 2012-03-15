\name{GQdk}
\alias{GQdk}
\title{Sparse Gaussian Quadrature grid}
\usage{
  GQdk(d = 1L, k = 1L)
}
\arguments{
  \item{d}{integer scalar - the dimension of the function
  to be integrated with respect to the standard
  \code{d}-dimensional Gaussian density}

  \item{k}{integer scalar - the order of the grid.  A grid
  of order \code{k} provides an exact result for a
  polynomial of total order of \code{2k - 1} or less
  multiplied by the}
}
\value{
  a matrix with \code{d + 1} columns.  The first column is
  the weights and the remaining \code{d} columns are the
  node coordinates.
}
\description{
  Generate the sparse multidimensional Gaussian quadrature
  grids
}
\note{
  The number of nodes gets very large very quickly with
  increasing \code{d} and \code{k}.  See the charts at
  \url{http://www.sparse-grids.de}.
}
\examples{
GQdk(2,5)
}

