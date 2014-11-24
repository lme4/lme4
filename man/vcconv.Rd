\name{vcconv}
\alias{vcconv}
\alias{mlist2vec}
\title{Convert between representations of (co-)variance structures}
\description{
  Convert between representations of (co-)variance structures.
}
\usage{
  mlist2vec(L)
}
\arguments{
  \item{L}{List of symmetric, upper-triangular, or lower-triangular
    square matrices.}
}
\value{
}
\details{
  \itemize{
    \item{\code{mlist2vec}}{Convert list of matrices to concatenated
      vector of lower triangles with an attribute that gives the
      dimension of each matrix in the original list.  This attribute may
      be used to reconstruct the matrices.}
  }
}
