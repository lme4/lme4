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
  The value differs based on the family.
  \itemize{
    \item For Gaussian fitted models (LMMs), this is the \dQuote{residual 
    standard deviation}.
    \item For Gamma fitted models (GLMMs), it represents the square root 
    of the dispersion parameter.
    \item For Binomial and Poisson fitted models, the value is always reported 
    as \code{1}, since their variance is determined entirely by the mean and 
    there is no separate scale parameter.
    %\item Less interpretably, the square root of the residual deviance per 
    %degree of freedom in more general models.
  }
}
\details{
  Package \CRANpkg{lme4} provides methods for mixed-effects models of class
  \code{\linkS4class{merMod}} and lists of linear models, 
  \code{\linkS4class{lmList4}}.
}
\examples{
methods(sigma)# from R 3.3.0 on, shows methods from pkgs 'stats' *and* 'lme4'
}
