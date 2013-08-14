\name{pirls}
\alias{pirls}
\alias{pirls1}
\title{Penalized iteratively reweighted least squares}
\usage{
  pirls(theta, beta, u, mu, eta, glmod, y,
    family = binomial, weights = rep(1, length(y)),
    offset = rep(0, length(y)), tol = 10^-6, npirls = 30,
    nAGQ = 0)

  pirls1(thetabeta, ...)
}
\arguments{
  \item{theta}{covariance parameters}

  \item{beta}{initial beta}

  \item{u}{initial u}

  \item{mu}{fitted values}

  \item{eta}{linear predictor}

  \item{glmod}{output of \code{glFormula}}

  \item{y}{response}

  \item{family}{family}

  \item{weights}{prior weights}

  \item{offset}{offset}

  \item{tol}{convergence tolerance}

  \item{npirls}{maximum number of iterations}

  \item{nAGQ}{either 0 (PIRLS for \code{u} and \code{beta})
  or 1 (\code{u} only).  currently no quadature is
  available}

  \item{thetabeta}{\code{c(theta,beta)} (\code{pirls1}
  only)}

  \item{...}{Arguments to pass to \code{pirls}}
}
\value{
  the laplace approximated deviance with a bunch of
  attributes: \itemize{ \item[pdev] penalized deviance at
  convergence \item[beta] fixed effects at convergence
  (\code{pirls0} only) \item[u] spherized random effects at
  convergence }
}
\description{
  A pure \code{R} implementation of Doug Bates' penalized
  iteratively reweighted least squares (PIRLS) algorithm
  for computing generalized linear mixed model deviances.
  The purpose is to clarify how PIRLS works without having
  to read through C++ code, and as a sandbox for trying out
  modified versions of PIRLS.
}
\details{
  \code{pirls1} is a convenience function for optimizing
  \code{pirls} under \code{nAGQ = 1}. In particular, it
  wraps \code{theta} and \code{beta} into a single argument
  \code{thetabeta}.
}

