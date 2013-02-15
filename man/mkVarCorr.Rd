\name{mkVarCorr}
\alias{mkVarCorr}
\title{Make variance and correlation matrices from \code{theta}}
\usage{
  mkVarCorr(sc, cnms, nc, theta, nms)
}
\arguments{
  \item{sc}{scale factor (residual standard deviation)}

  \item{cnms}{component names}

  \item{nc}{numeric vector: number of terms in each RE
  component}

  \item{theta}{theta vector (lower-triangle of Cholesky
  factors)}

  \item{nms}{component names (FIXME: nms/cnms redundant:
  nms=names(cnms)?)}
}
\value{
  A matrix
}
\description{
  Make variance and correlation matrices from \code{theta}
}
\seealso{
  \code{\link{VarCorr}}
}

