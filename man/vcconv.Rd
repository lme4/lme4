\name{vcconv}
\alias{vcconv}
\alias{mlist2vec}
\alias{vec2mlist}
\alias{vec2STlist}
\alias{sdcor2cov}
\alias{cov2sdcor}
\alias{Vv_to_Cv}
\alias{Sv_to_Cv}
\alias{Cv_to_Vv}
\alias{Cv_to_Sv}
\title{Convert between representations of (co-)variance structures}
\description{
  Convert between representations of (co-)variance structures
  (EXPERIMENTAL).  See source code for details.
}
\usage{
  mlist2vec(L)
  vec2mlist(v, n = NULL, symm = TRUE)
  vec2STlist(v, n = NULL)
  sdcor2cov(m)
  cov2sdcor(V)
  Vv_to_Cv(v, n = NULL, s = 1)
  Sv_to_Cv(v, n = NULL, s = 1)
  Cv_to_Vv(v, n = NULL, s = 1)
  Cv_to_Sv(v, n = NULL, s = 1)
}
\arguments{
  \item{L}{List of symmetric, upper-triangular, or lower-triangular
    square matrices.}
  \item{v}{Concatenated vector containing the elements of the
    lower-triangle (including the diagonal) of a symmetric or triangular
    matrix.}
  \item{n}{Number of rows (and columns) of the resulting matrix.}
  \item{symm}{Return symmetric matrix if \code{TRUE} or lower-triangular
    if \code{FALSE}.}
  \item{m}{Standard deviation-correlation matrix.}
  \item{V}{Covariance matrix.}
  \item{s}{Scale parameter.}
}
\value{
  (Co-)variance structure
}
\details{
  \itemize{
    \item{\code{mlist2vec}}{Convert list of matrices to concatenated
      vector of lower triangles with an attribute that gives the
      dimension of each matrix in the original list. This attribute may
      be used to reconstruct the matrices.  Returns a concatenation of
      the elements in one triangle of each matrix. An attribute
      \code{"clen"} gives the dimension of each matrix.}
    \item{\code{vec2mlist}}{Convert concatenated vector to list of matrices (lower
      triangle or symmetric). These matrices could represent Cholesky
      factors, covariance matrices, or correlation matrices (with
      standard deviations on the diagonal).}
    \item{\code{vec2STlist}}{Convert concatenated vector to list of ST
      matrices.}
    \item{\code{sdcor2cov}}{Standard deviation-correlation matrix to
      covariance matrix convert 'sdcor' format (std dev on diagonal,
      cor on off-diag) to and from variance-covariance matrix.}
    \item{\code{cov2sdcor}}{Covariance matrix to standard
      deviation-correlation matrix (i.e. standard deviations on the
      diagonal and correlations off the diagonal).}
    \item{\code{Vv_to_Cv}}{Variance-covariance to relative covariance
      factor. Returns a vector of elements from the lower triangle of a
      relative covariance factor.}
    \item{\code{Sv_to_Cv}}{Standard-deviation-correlation to relative covariance
      factor. Returns a vector of elements from the lower triangle of a
      relative covariance factor.}
    \item{\code{Cv_to_Vv}}{Relative covariance factor to
      variance-covariance. From unscaled Cholesky vector to (possibly
      scaled) variance-covariance vector. Returns a vector of elements
      from the lower triangle of a variance-covariance matrix.}
    \item{\code{Cv_to_Sv}}{Relative covariance factor to
      standard-deviation-correlation.  From unscaled Chol to sd-cor
      vector. Returns a vector of elements from the lower triangle of a
      standard-deviation-correlation matrix.}
  }
}
\examples{
vec2mlist(1:6)
mlist2vec(vec2mlist(1:6)) # approximate inverse
}
