\name{getL}
\title{Extract the  L  Matrix from a *Mer object}
\alias{getL}
\alias{getL-methods}
\alias{getL,merMod-method}
\alias{getL,reModule-method}
\description{
  Extract the Cholesky factor \eqn{L} from a \code{*mer} fitted object.
  A simple utility returning \eqn{L} independently from the inner
  detailed structure of \code{x}.
}
\usage{
getL(x)
}
\arguments{
  \item{x}{a fitted model, result of \code{\link{lmer}()},
    \code{\link{glmer}()}, or similar.}
}
\section{Methods}{
  \describe{
    \item{\code{signature(x = "lmerMod")}}{ }%% FIXME
    \item{\code{signature(x = "mer")}}{ }
    \item{\code{signature(x = "merenv")}}{ }
    \item{\code{signature(x = "reModule")}}{ }
  }
}
\details{
  The Cholesky factor \eqn{L} is a lower triangular
  \eqn{q \times q}{q x q} matrix, defined by
  \deqn{L L' = \Lambda' Z'Z \Lambda + I_q}{L L' = Lambda' Z'Z Lambda + I_q}, 
  where \eqn{I_q} is the identity and
  \eqn{\Lambda = \Lambda_\theta = \Lambda(\theta)}{Lambda = Lambda(\theta)}
  is the relative covariance factor (typically lower triangular as
  well), generating the \eqn{q \times q}{q x q} variance-covariance
  matrix \eqn{\Sigma_\theta}{Sigma(theta)} as
  \deqn{\Sigma_\theta = \sigma^2 \Lambda_\theta{\Lambda_{\theta}}^T.}{Sigma = sigma^2 Lambda Lambda^T.}
}
\value{
  the Cholesky factor \eqn{L}, an object extending class
  \code{\linkS4class{CHMfactor}},
  (currently typically of class \code{\linkS4class{CHMsimpl}}.
}
%\references{Douglas Bates (2010).  lme4: Mixed Models with R
%}
% \seealso{
% }
\examples{
## (see  help(Penicillin)):
fm1 <- lmer(diameter ~ (1|plate) + (1|sample), Penicillin)
image(getL(fm1),
      sub = "Structure of random effects interaction in penicillin model")
}
\keyword{methods}
