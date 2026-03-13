\name{sigma}
\title{Extract Residual Standard Deviation 'Sigma'}
\alias{sigma}
\alias{sigma.merMod}
\description{
  Extract the estimated standard deviation of the errors, the
  \dQuote{residual standard deviation} (also misnamed the
  \dQuote{residual standard error}), from a fitted model of class
  \code{\linkS4class{merMod}}.
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
  The value differs based on the family (see details for an explanation.)
  \itemize{
    \item For LMMs (Gaussian family), this is the residual standard deviation
    (ML or REML estimate).
    \item For GLMMs with a scale parameter (e.g. Gamma, Inverse Gaussian),
    it is \code{sqrt(deviance(object) / df.residual(object))}, consistent
    with the GLM convention.
    \item For GLMMs without a separate scale parameter (Binomial, Poisson),
    the value is always \code{1}.
  }
}
\details{
  In general the dispersion parameter (which we call \code{sigma}) is the 
  square root of the constant multiplier in the 'variance function' provided 
  by the \code{family} functions.
  \itemize{
    \item For LMMs (Gaussian family), \code{gaussian()$variance = rep.int(1, length(mu))}
    and the variance is \eqn{\sigma^{2}}. The ML estimate is
    \eqn{\hat\sigma = \sqrt{\text{pwrss}/n}}, and the REML estimate uses
    \eqn{n - p} in the denominator.
    \item For GLMMs with a scale parameter, \code{sigma} is estimated as
    \eqn{\sqrt{D / \nu}} where \eqn{D} is the sum of squared deviance
    residuals (\code{deviance(object)}) and \eqn{\nu} is the residual degrees
    of freedom (\code{df.residual(object)}). This is consistent with how
    \code{stats::sigma} works for \code{glm} objects.
    \item For Gamma fits, \code{Gamma()$variance = mu^2} where 
    \code{mu = scale*shape}. The known variance is 
    \eqn{\frac{\mu^{2}}{\code{shape}}}.
    Thus, the constant multiplier here is \eqn{\frac{1}{\sqrt{\code{shape}}}},
    which \code{sigma} estimates as \eqn{\sqrt{D / \nu}}.
    \item Similarly, for Inverse Gaussian fits, we have
    \code{inverse.gaussian()$variance = mu^{3}}, with known variance
    \eqn{\frac{\mu^{3}}{\lambda}}; the constant multiplier here is 
    \eqn{\frac{1}{\sqrt{\lambda}}} (\eqn{\lambda} is the shape parameter),
    estimated as \eqn{\sqrt{D / \nu}}.
  }
}
\seealso{
  Package \CRANpkg{lme4} provides methods for mixed-effects models of class
  \code{\linkS4class{merMod}} and lists of linear models, 
  \code{\linkS4class{lmList4}}.
}
\examples{
methods(sigma)# from R 3.3.0 on, shows methods from pkgs 'stats' *and* 'lme4'
}
