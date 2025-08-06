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
  The value differs based on the family (see details above for an explanation.)
  \itemize{
    \item For Gaussian fitted models, this is the residual 
    standard deviation.
    \item For Gamma and Inverse Gaussian fitted models, it represents 
    the square root of the inverse of the shape parameter.
    \item For Binomial and Poisson fitted models, the value is always reported 
    as \code{1}, since their variance is determined entirely by the mean and 
    there is no separate scale parameter.
    %\item Less interpretably, the square root of the residual deviance per 
    %degree of freedom in more general models.
  }
}
\details{
  In general the dispersion parameter (which we call \code{sigma}) is the 
  square root of the constant multiplier in the 'variance function' provided 
  by the \code{family} functions.
  \itemize{
    \item For Gaussian fits, \code{gaussian()$variance = rep.int(1, length(mu))} 
    and it's known that the variance is \eqn{\sigma^{2}}. Hence, the constant
    multiplier of the variance function in this case is \eqn{\sigma^{2}}, so 
    we say \code{sigma} is \eqn{\sqrt{\sigma}}.
    \item For Gamma fits, \code{Gamma()$variance = mu^2} where 
    \code{mu = scale*shape}. The known variance is 
    \eqn{\frac{\mu^{2}}{\code{shape}}}.
    Thus, the constant multiplier here is \eqn{\frac{1}{\sqrt{\code{shape}}}}.
    \item Similarly, for Inverse Gaussian fits, we have
    \code{inverse.gaussian()$variance = mu^{3}}, with known variance
    \eqn{\frac{\mu^{3}}{\lambda}}; similarly, the constant multiplier here is 
    \eqn{\frac{1}{\sqrt{\lambda}}}. \eqn{\lambda} is referred to as the shape
    parameter.
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
