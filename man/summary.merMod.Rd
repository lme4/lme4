\name{summary}
\alias{summary}
\alias{summary.merMod}
\alias{show.summary.merMod}
\alias{print.summary.merMod}
\title{
  Summary for a [ng]lmer Fit
}
\description{
  \code{summary.merMod} Computes and returns a list of summary statistics of a 
  fitted mixed-effects model of class \code{\linkS4class{merMod}}, as produced 
  by [ng]lmer fits. The summary includes estimates for the fixed and random 
  effects. \code{print.summary} controls the output for the summary method.
}
\usage{
\S3method{summary}{merMod}(object,
       correlation = (p <= getOption("lme4.summary.cor.max")),
       use.hessian = NULL, \dots)
\S3method{print}{summary.merMod}(x, digits = max(3, getOption("digits") - 3),
      correlation = NULL, symbolic.cor = FALSE,
      signif.stars = getOption("show.signif.stars"),
      ranef.comp = c("Variance", "Std.Dev."),
      ranef.corr = any(ranef.comp == "Std.Dev."), show.resids = TRUE, \dots)
}
\arguments{
  \item{object}{a fitted [ng]lmer model.}
  \item{correlation}{(logical) indicates whether the correlation matrix
    should be computed and stored along with the covariance.
    For \code{print.summary.merMod}, indicates whether the correlation
    matrix of the fixed-effects parameters should be printed.  In the
    latter case, when \code{NULL} (the default), the correlation matrix
    is printed when it has been computed by \code{summary(.)}, and when
    \eqn{p <= 12}, and the cutoff \code{12} may be modified by
    \code{options(lme4.summary.cor.max = <n>)}}
  \item{use.hessian}{(logical) indicates whether to use the
    finite-difference Hessian of the deviance function to compute
    standard errors of the fixed effects; see \code{\link{vcov.merMod}}
    for details}
  \item{x}{an \R object of class \code{summary.merMod}.}
  \item{digits}{number of significant digits for printing.}
  \item{symbolic.cor}{should a symbolic encoding of the fixed-effects 
    correlation matrix be printed?  If so, the \code{\link{symnum}} function is 
    used.}
  \item{signif.stars}{(logical) should significance stars be used?}
  \item{ranef.comp}{character vector of length one or two, indicating
    if random-effects parameters should be reported on the variance and/or
    standard deviation scale.}
  \item{ranef.corr}{(logical) print correlations (rather than covariances) of
    random effects?}
  \item{show.resids}{should the quantiles of the scaled residuals be printed?}
  \item{\dots}{potentially further arguments passed from other methods.}
}
\seealso{
  \itemize{
    \item\code{\link[base]{summary}} for the original function in base \code{R},
    \item\code{\link[lme4]{lmer}}, and \code{\link[lme4]{glmer}} for model
      fitting,
    \item\code{\link[lme4]{sigma.merMod}} documents the result obtained in
      \code{"Residual Std.Dev"}; the name "Residual Std.Dev" is mostly 
      applicable for \code{\link[lme4]{lmer}} fits, but the meaning differs for
      \code{\link[lme4]{glmer}} fits.
  }
}
\value{
  An S3 object of class \code{"summary.merMod"}, which is a list containing 
  fixed and random effect estimates along with other standard model summary 
  information.
}
\examples{
fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
smry <- summary(fm1)
# Obtaining the variance-covariance matrix of the fixed effects
smry$vcov
# Obtaining the correlation matrix of fixed effects
smry$vcov@factors$correlation
}