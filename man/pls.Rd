\name{pls}
\alias{pls}
\alias{pls.list}
\alias{pls.matrix}
\title{Create linear mixed model deviance function}
\usage{
  pls(obj, theta, y, ...)

  \method{pls}{matrix} (X, y, Zt, Lambdat, thfun,
    weights = rep(1, n), offset = numeric(n), REML = TRUE)

  \method{pls}{list} (obj, y, weights = rep(1, length(y)),
    offset = rep(0, length(y)), REML = TRUE)
}
\arguments{
  \item{obj}{output of \code{lFormula} or a model matrix}

  \item{y}{response}

  \item{weights}{prior weights}

  \item{offset}{offset}

  \item{REML}{calculate REML deviance?}
}
\value{
  a function that evaluates the deviance or REML criterion
}
\description{
  A pure \code{R} implementation of the penalized least
  squares (PLS) approach for computing linear mixed model
  deviances. The purpose is to clarify how PLS works
  without having to read through C++ code, and as a sandbox
  for trying out modifications to PLS.
}

