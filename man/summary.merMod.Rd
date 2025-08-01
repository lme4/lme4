\name{summary}
\alias{summary}
\alias{summary.merMod}
\title{
  Summary for a [ng]lmer Fit
}
\description{
  Summarizes the results of a fitted mixed-effects model of class 
  \code{\linkS4class{merMod}}, as produced by [ng]lmer fits. The summary 
  includes estimates for the fixed and random effects. 
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