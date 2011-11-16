\name{lmResp}
\alias{lmResp}
\alias{lmerResp}
\alias{glmResp}
\alias{nlsResp}
\title{Generator objects for the response classes}
\section{Methods}{
  \describe{
    \item{\code{new(y=y)}:}{Create a new
      \code{\linkS4class{lmResp}} or \code{\linkS4class{lmerResp}} object}
    \item{\code{new(family=family, y=y)}:}{Create a new
      \code{\linkS4class{glmResp}} object}
    \item{\code{new(y=y, nlmod=nlmod, nlenv=nlenv, pnames=pnames, N=N)}:}{Create a new
      \code{\linkS4class{nlsResp}} object}
  }
}
\arguments{
  \item{y}{the numeric response vector}
  \item{family}{a \code{\link{family}} object}
  \item{nlmod}{the nonlinear model function}
  \item{nlenv}{an environment holding data objects for evaluation of
    \code{nlmod}}
  \item{pnames}{a character vector of parameter names}
  \item{N}{an integer scalar - the length of the linear predictor}
}
\note{Arguments to the \code{new} methods must be named arguments.}
\description{The generator objects for the \code{\linkS4class{lmResp}},
  \code{\linkS4class{lmerResp}}, \code{\linkS4class{glmResp}} and
  \code{\linkS4class{nlsResp}} reference classes. Such objects are
  primarily used through their \code{new} methods.}
\seealso{
  \code{\linkS4class{lmResp}},
  \code{\linkS4class{lmerResp}},
  \code{\linkS4class{glmResp}},
  \code{\linkS4class{nlsResp}}
}
\keyword{classes}
