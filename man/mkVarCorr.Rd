\name{mkVarCorr}
\title{Make Variance and Correlation Components of Covariance}
\alias{mkVarCorr}
\usage{
mkVarCorr(sc, cnms, nc, theta, nms, reCovs = NULL,
          full_cor = NULL, is_lmm = NULL)
}
\description{
  Computes from \code{reCovs} the variance and correlation components of
  random effect covariance matrices.
}
\arguments{
  \item{sc}{a numeric vector of length 1 giving the residual standard
    deviation for LMMs.}
  \item{cnms}{a list of character vectors giving component dimnames.}
  \item{nc}{an integer vector giving component dimensions (number of
    columns).}
  \item{theta}{a numeric vector giving covariance parameters
    (lower triangular entries of relative Cholesky factors in
    column-major order), used only to determine \code{reCovs} when
    that is missing.}
  \item{nms}{a character vector giving component names.}
  \item{reCovs}{a list of \code{\linkS4class{Covariance}} objects.}
  \item{full_cor}{a logical indicating if correlation matrices are
    computed.}
  \item{is_lmm}{a logical indicating if the model is an LMM as opposed
    to a GLMM, determining if \code{sc} is used.}
}
\value{
  A list of matrices.
}
\seealso{\code{\link{VarCorr}}.}
