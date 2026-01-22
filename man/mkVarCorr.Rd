\name{mkVarCorr}
\title{Make Variance and Correlation Matrices from \code{theta}}
\alias{mkVarCorr}
\usage{
mkVarCorr(sc, cnms, nc, theta, nms, reCovs = NULL, full_cor = NULL, is_lmm = NULL)
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
  \item{reCovs}{a list of \code{Covariance} objects.}
  \item{full_cor}{logical; specifies whether the full correlation matrix should be 
    produced.}
  \item{is_lmm}{logical; indicating whether the model is an LMM as GLMMs do 
    not include a residual scale parameter.}
  
}
\value{
  A \code{\link{matrix}}
}
\seealso{
  \code{\link{VarCorr}}
}

