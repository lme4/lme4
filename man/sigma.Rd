\name{sigma}
\alias{sigma}
\alias{sigma.merMod}
\title{Extract residual standard error}
\usage{
\S3method{sigma}{merMod}(object, \dots)
}
\arguments{
  \item{object}{a fitted model.}

  \item{\dots}{additional, optional arguments.  (None are
  used in the merMod method)}
}
\value{
  the residual standard error as a scalar
}
\description{
  Extract the residual standard error from a fitted model.
}
\details{
  Package \pkg{lme4} provides methods for mixed-effects models of class
  \code{\linkS4class{merMod}} and lists of LM's, \code{\linkS4class{lmList4}}.
}
\examples{
methods(sigma)# from R 3.3.0 on, shows methods from pkgs 'stats' *and* 'lme4'
}
