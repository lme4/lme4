\name{pls}
\alias{pls}
\alias{pls.list}
\alias{pls.matrix}
\title{Create linear mixed model deviance function}
\usage{
  pls(obj, y, ...)

  \method{pls}{matrix} (obj, y, Zt, Lambdat, thfun,
    weights = rep(1, n), offset = numeric(n), REML = TRUE,
    ...)

  \method{pls}{list} (obj, y, weights = rep(1, length(y)),
    offset = rep(0, length(y)), REML = TRUE, ...)
}
\arguments{
  \item{obj}{output of \code{lFormula} or a model matrix}

  \item{y}{response}

  \item{...}{Arguments to pass to other functions}

  \item{Zt}{transpose of the sparse model matrix for the
  random effects}

  \item{Lambdat}{upper triangular sparse Cholesky factor of
  the relative covariance matrix of the random effects}

  \item{thfun}{a function that takes a value of
  \code{theta} and produces the non-zero elements of
  \code{Lambdat}.  The structure of \code{Lambdat} cannot
  change, only the numerical values}

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
\examples{
library(lme4)
library(nloptwrap)
lmod <- lFormula(Reaction ~ Days + (Days|Subject), sleepstudy)
devf <- pls(lmod,sleepstudy$Reaction)
devf(c(1,0,1))             # starting value
bobyqa(c(1, 0, 1), devf, lower=c(0,-Inf,0))[c("par","value")]
}
\keyword{models}

