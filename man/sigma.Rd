\name{sigma}
\title{Extract Residual Standard Deviation 'Sigma'}
\alias{sigma}
\alias{sigma.merMod}
\description{
  Extract the estimated standard deviation of the errors, the
  \dQuote{residual standard deviation} (also misnamed the
  \dQuote{residual standard error}), from a fitted model.
}
\usage{
\S3method{sigma}{merMod}(object, \dots)
}
\arguments{
  \item{object}{a fitted model.}

  \item{\dots}{additional, optional arguments, passed from or to
    methods. (None currently in our two methods.)}
}
\value{
  Typically a number, the estimated standard deviation of the
  errors (\dQuote{residual standard deviation}) for Gaussian
  models, and - less interpretably - the square root of the residual
  deviance per degree of freedom in more general models.
}
\details{
  Package \pkg{lme4} provides methods for mixed-effects models of class
  \code{\linkS4class{merMod}} and lists of linear models, \code{\linkS4class{lmList4}}.
}
\examples{
methods(sigma)# from R 3.3.0 on, shows methods from pkgs 'stats' *and* 'lme4'
}
