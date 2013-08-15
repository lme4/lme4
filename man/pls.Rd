\name{pls}
\alias{pls}
\alias{pls.list}
\alias{pls.matrix}
\title{Compute linear mixed model deviance}
\usage{
  pls(obj, theta, y, ...)

  \method{pls}{matrix} (obj, theta, y, Zt, Lambdat, Lind,
    weights = rep(1, length(y)),
    offset = rep(0, length(y)), REML = TRUE)

  \method{pls}{list} (obj, theta, y,
    weights = rep(1, length(y)),
    offset = rep(0, length(y)), REML = TRUE)
}
\arguments{
  \item{obj}{output of \code{lFormula} or a model matrix}

  \item{theta}{covariance parameters}

  \item{y}{response}

  \item{weights}{prior weights}

  \item{offset}{offset}

  \item{REML}{calculate REML deviance?}
}
\value{
  deviance with a bunch of attributes: \itemize{
  \item[devML] ML deviance \item[devREML] REML deviance
  \item[beta] fixed effects at convergence \item[u]
  conditional mode at convergence \item[sigmaML] ML
  estimate of the residual standard deviation
  \item[sigmaREML] REML estimate of the residual standard
  deviation }
}
\description{
  A pure \code{R} implementation of Doug Bates' penalized
  least squares (PLS) approach for computing linear mixed
  model deviances. The purpose is to clarify how PLS works
  without having to read through C++ code, and as a sandbox
  for trying out modifications to PLS.
}

