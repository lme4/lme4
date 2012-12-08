\name{RSCupdate}
\alias{RSCupdate}
\title{Update for the penalized least squares problem}
\usage{
  RSCupdate(rv, xv, theta, m)
}
\arguments{
  \item{rv}{the matrix of row indices for the regular
  sparse column representation of ZtXt}

  \item{xv}{the non-zero values in ZtXt}

  \item{theta}{vector of relative standard deviations}

  \item{m}{the sparse matrix to be updated.  Must have the
  correct pattern.}
}
\description{
  Update for the penalized least squares problem
}
\examples{
rv <- rbind(as.integer(Dyestuff$Batch) - 1L, 6L)
xv <- matrix(1, nrow=nrow(rv), ncol=ncol(rv))
A <- sparseMatrix(i=as.integer(c(1,1,2,2,3,3,4,4,5,5,6,6,7)),
                  j=as.integer(c(1,7,2,7,3,7,4,7,5,7,6,7,7)),
                  x=rep.int(1.,13L), symmetric=TRUE,check=FALSE)
RSCupdate(rv, xv, 0.831, A)
A
}

