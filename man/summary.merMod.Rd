\name{summary.merMod}
\alias{summary.merMod}
\title{
  Summary for a [ng]lmer Fit
}
\description{
  Summarizes the results of a fitted mixed-effects model of class 
  \code{merMod}, as produced by [ng]lmer fits. The summary includes estimates
  for the fixed and random effects. 
}
\usage{
\method{summary}{merMod}(object,
       correlation = (p <= getOption("lme4.summary.cor.max")),
       use.hessian = NULL, \dots)
}
\arguments{
  \item{object}{a fitted [ng]lmer model.}
  \item{correlation}{(logical) if \code{TRUE}, the correlation matrix of the 
    fixed effects is included in the output.}
  \item{use.hessian}{(logical) if \code{TRUE}, uses the Hessian to compute the 
  variance-covariance matrix of the fixed effects.}
}
}
\details{
  \itemize{
    \item For Gamma GLMMs, the \code{"Residual Std.Dev."} corresponds to the 
    square root of the dispersion parameter rather than a residual standard 
    deviation.
    \item TODO: add info for other GLMMs....
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