\name{isREML}
\alias{isGLMM}
\alias{isLMM}
\alias{isNLMM}
\alias{isREML}
\title{Check characteristics of models}
\usage{
  isREML(x, ...)

  isLMM(x, ...)

  isNLMM(x, ...)

  isGLMM(x, ...)
}
\arguments{
  \item{x}{a fitted model.}

  \item{...}{additional, optional arguments.  (None are
  used in the merMod methods)}
}
\value{
  a logical value
}
\description{
  Check characteristics of models: whether a model fit
  corresponds to a linear (LMM), generalized linear (GLMM),
  or nonlinear (NLMM) mixed model, and whether a linear
  mixed model has been fitted by REML or not
  (\code{isREML(x)} is always \code{FALSE} for GLMMs and
  NLMMs).
}
\details{
  This is a generic function.  At present the only methods
  are for mixed-effects models of class
  \code{\linkS4class{merMod}}.
}
\seealso{
  getME
}

