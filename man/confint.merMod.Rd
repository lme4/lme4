\name{confint.merMod}
\alias{confint.merMod}
\title{Compute confidence intervals on the parameters of an lme4 fit}
\usage{
  confint.merMod(object, parm, level = 0.95,
    method = c("profile", "Wald", "boot"), zeta,
    nsim = 500, boot.type = "perc", quiet = FALSE, ...)
}
\arguments{
  \item{object}{a fitted [ng]lmer model}

  \item{parm}{parameters (specified by integer position)}

  \item{level}{confidence level}

  \item{method}{for computing confidence intervals}

  \item{zeta}{likelihood cutoff (if not specified, computed
  from \code{level}: "profile" only)}

  \item{nsim}{number of simulations for parametric
  bootstrap intervals}

  \item{boot.type}{bootstrap confidence interval type}

  \item{quiet}{(logical) suppress messages about
  computationally intensive profiling?}

  \item{\dots}{additional parameters to be passed to
  \code{\link{profile.merMod}}}
}
\value{
  a numeric table of confidence intervals
}
\description{
  Compute confidence intervals on the parameters of an lme4
  fit
}
\details{
  Depending on the method specified, this function will
  compute confidence intervals by ("profile") computing a
  likelihood profile and finding the appropriate cutoffs
  based on the likelihood ratio test; ("Wald") approximate
  the confidence intervals (of fixed-effect parameters
  only) based on the estimated local curvature of the
  likelihood surface; ("boot") perform parametric
  bootstrapping with confidence intervals computed from the
  bootstrap distribution according to \code{boot.type} (see
  \code{\link{boot.ci}})
}

