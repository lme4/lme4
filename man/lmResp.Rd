\name{lmResp}
\alias{glmResp}
\alias{lmerResp}
\alias{lmResp}
\alias{nlsResp}
\title{Generator objects for the response classes}
\arguments{
  \item{y}{the numeric response vector}

  \item{family}{a \code{\link{family}} object}

  \item{nlmod}{the nonlinear model function}

  \item{nlenv}{an environment holding data objects for
  evaluation of \code{nlmod}}

  \item{pnames}{a character vector of parameter names}

  \item{gam}{a numeric vector - the initial linear
  predictor}
}
\description{
  The generator objects for the
  \code{\linkS4class{lmResp}},
  \code{\linkS4class{lmerResp}},
  \code{\linkS4class{glmResp}} and
  \code{\linkS4class{nlsResp}} reference classes. Such
  objects are primarily used through their \code{new}
  methods.
}
\note{
  Arguments to the \code{new} methods must be named
  arguments.
}
\section{Methods}{
  \describe{ \item{\code{new(y=y)}:}{Create a new
  \code{\linkS4class{lmResp}} or
  \code{\linkS4class{lmerResp}} object.}
  \item{\code{new(family=family, y=y)}:}{Create a new
  \code{\linkS4class{glmResp}} object.}
  \item{\code{new(y=y, nlmod=nlmod, nlenv=nlenv,
  pnames=pnames, gam=gam)}:}{Create a new
  \code{\linkS4class{nlsResp}} object.} }
}
\seealso{
  \code{\linkS4class{lmResp}},
  \code{\linkS4class{lmerResp}},
  \code{\linkS4class{glmResp}},
  \code{\linkS4class{nlsResp}}
}
\keyword{classes}

