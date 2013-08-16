\name{mkReTrm}
\alias{mkReTrm}
\title{Generate a random effect design from a standard '|' term}
\usage{
  mkReTrm(bar, fr)
}
\arguments{
  \item{bar}{a random effect specification as returned by
  \code{\link{findbars}}.}

  \item{fr}{a \code{model.frame}}
}
\value{
  a list with components \item{ff}{the grouping factor}
  \item{Zt}{transpose of the sparse model matrix for the
  random effect} \item{Lambdat}{transpose of the sparse
  relative covariance factor} \item{theta}{initial values
  of the covariance parameters} \item{Lind}{index vector:
  which entries in \code{theta} go into which entries in
  the \code{x} slot of \code{Lambdat} (included for
  backwards compatibility)} \item{thfun}{a function mapping
  theta to the \code{"x"} slot of \code{Lambdat} -- in this
  case, simply \code{theta[Lind]}} \item{lower}{lower
  bounds on the covariance parameters} \item{upper}{upper
  bounds on the covariance parameters} \item{cnms}{a list
  of column names of the random effects according to the
  grouping factors} \item{special}{\code{FALSE}}
  \item{nlambda}{the number of non-zero entries in
  \code{Lambdat}} \item{ntheta}{\code{length(theta)},
  number of (co-)variance parameters} \item{nc}{the number
  of ranefs per level} \item{nb}{\code{nrow(Zt)}, the
  number of random effects} \item{nl}{number of levels of
  this grouping factor}
}
\description{
  Generates the partial design matrix, the partial Cholesky
  factor for an unstrucured covariance of the random
  effects, a function mapping the variance-covariance
  parameters to the entries in the partial Cholesky factor
  and additional information required to construct the full
  models for one random effect term in the .
}

