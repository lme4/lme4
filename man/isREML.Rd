\name{isREML}
\alias{isREML}
\title{Check for a REML fit}
\usage{
  isREML(x, ...)
}
\arguments{
  \item{x}{a fitted model.}

  \item{...}{additional, optional arguments.  (None are
  used in the merMod method)}
}
\value{
  \code{TRUE} if \code{x} has been fit by REML, otherwise
  FALSE
}
\description{
  Check if a model has been fit according to the REML
  criterion
}
\details{
  This is a generic function.  At present the only methods
  are for mixed-effects models of class
  \code{\linkS4class{merMod}}.
}

