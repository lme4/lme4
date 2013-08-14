\name{pls}
\alias{pls}
\title{Compute linear mixed model deviance}
\usage{
  pls(theta, lmod, y, weights = rep(1, length(y)),
    offset = rep(0, length(y)), REML = TRUE)
}
\arguments{
  \item{theta}{covariance parameters}

  \item{lmod}{output of \code{lFormula}}

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

