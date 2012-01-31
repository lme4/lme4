\docType{data}
\name{GQN}
\alias{GQN}
\title{Sparse Gauss-Hermite quadrature grids}
\format{A list of lists.}
\description{
  \code{GQN} contains the non-redundant quadrature nodes
  and weights for integration of a scalar function of a
  \code{d}-dimensional argument with respect to the density
  function of the \code{d}-dimensional Gaussian density
  function.  These are stored in a list of lists.  The
  outer list is indexed by the dimension, \code{d}, in the
  range of 1 to 20.  The inner list is indexed by \code{k},
  the order of the quadrature.
}
\note{
  These are only the non-redundant nodes.  To regenerate
  the whole array of nodes, all possible permutations of
  axes and all possible combinations of \eqn{\pm 1}{+/- 1}
  must be applied to the axes. The function
  \code{\link{GQdk}} reproduces the entire array of nodes.
}
\examples{
GQN[[3]][[5]]
}
\seealso{
  \code{\link{GQdk}}
}
\keyword{datasets}

