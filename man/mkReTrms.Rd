\name{mkReTrms}
\alias{mkReTrms}
\title{Create Z, Lambda, Lind, etc.}
\usage{
  mkReTrms(bars, fr)
}
\arguments{
  \item{bars}{a list of parsed random-effects terms}

  \item{fr}{a model frame in which to evaluate these terms}
}
\value{
  a list with components \item{Zt}{transpose of the sparse
  model matrix for the random effects}
  \item{Lambdat}{transpose of the sparse relative
  covariance factor} \item{Lind}{an integer vector of
  indices determining the mapping of the elements of the
  \code{theta} to the \code{"x"} slot of \code{Lambdat}}
  \item{theta}{initial values of the covariance parameters}
  \item{lower}{lower bounds on the covariance parameters}
  \item{flist}{list of grouping factors used in the
  random-effects terms} \item{cnms}{a list of column names
  of the random effects according to the grouping factors}
}
\description{
  From the result of \code{\link{findbars}} applied to a
  model formula and and the evaluation frame, create the
  model matrix, etc. associated with random-effects terms.
  See the description of the returned value for a detailed
  list.
}
\seealso{
  Other utilities: \code{\link{findbars}},
  \code{\link{mkRespMod}}, \code{\link{nlformula}},
  \code{\link{nobars}}, \code{\link{subbars}}
}

