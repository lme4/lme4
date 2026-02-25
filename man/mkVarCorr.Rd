\name{mkVarCorr}
\title{Make Variance and Correlation Components of Covariance}
\alias{mkVarCorr}
\usage{
mkVarCorr(sc, cnms, nc = lengths(cnms, use.names = FALSE),
          theta, nms = names(cnms), reCovs = NULL)
}
\description{
  Computes from \code{reCovs} the variance and correlation components of
  random effect covariance matrices.
}
\arguments{
  \item{sc}{a numeric vector of length 1 giving the residual standard
    deviation for LMMs.  Set to \code{NULL} for GLMMs.}
  \item{cnms}{a list of character vectors giving component dimnames,
    used only to determine \code{nc} and \code{nms} if those are
    missing.}
  \item{nc}{an integer vector giving component dimensions (number of
    columns).}
  \item{theta}{a numeric vector giving covariance parameters
    (lower triangular entries of relative Cholesky factors in
    column-major order), used only to determine \code{reCovs} when
    that is missing.}
  \item{nms}{a character vector giving component names.}
  \item{reCovs}{a list of \code{\linkS4class{Covariance}} objects.}
}
\value{
  A list of matrices.
}
\seealso{\code{\link{VarCorr}}.}
