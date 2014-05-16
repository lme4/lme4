\name{sigma}
\alias{sigma}
\alias{sigma.merMod}
\title{Extract residual standard error}
\usage{
  sigma(object, ...)
}
\arguments{
  \item{object}{a fitted model.}

  \item{...}{additional, optional arguments.  (None are
  used in the merMod method)}
}
\value{
  the residual standard error as a scalar
}
\description{
  Extract the residual standard error from a fitted model.
}
\details{
  This is a generic function.  At present the only methods
  are for mixed-effects models of class
  \code{\linkS4class{merMod}}.
}

