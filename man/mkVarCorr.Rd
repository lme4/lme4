\name{mkVarCorr}
\title{Make Variance and Correlation Matrices from \code{theta}}
\alias{mkVarCorr}
\usage{
mkVarCorr(sc, cnms, nc, theta, nms)
}
\description{
  Make variance and correlation matrices from \code{theta}
}
\arguments{
  \item{sc}{scale factor (residual standard deviation).}
  \item{cnms}{component names.}
  \item{nc}{numeric vector: number of terms in each RE component.}
  \item{theta}{theta vector (lower-triangle of Cholesky factors).}
  \item{nms}{component names (FIXME: nms/cnms redundant:
  nms=names(cnms)?)}
}
\value{
  A \code{\link{matrix}}
}
\seealso{
  \code{\link{VarCorr}}
}

