\name{confint.merMod}
\alias{confint.merMod}
\title{Compute confidence intervals on the parameters of an lme4 fit}
\usage{
\method{confint}{merMod}(object, parm, level = 0.95,
	method = c("profile", "Wald", "boot"), zeta,
	nsim = 500, boot.type = "perc", quiet = FALSE,
	oldNames = TRUE, ...)
}
\arguments{
  \item{object}{a fitted [ng]lmer model}
  \item{parm}{parameters (specified by integer position)}
  \item{level}{confidence level \eqn{< 1}, typically above 0.90.}
  \item{method}{a \code{\link{character}} string determining the method
    for computing the confidence intervals.}
  \item{zeta}{(for \code{method = "profile"} only:) likelihood cutoff
    (if not specified, as by default, computed from \code{level}).}
  \item{nsim}{number of simulations for parametric bootstrap intervals.}
  \item{boot.type}{bootstrap confidence interval type.}
  \item{quiet}{(logical) suppress messages about computationally intensive profiling?}
  \item{oldNames}{(logical) use old-style names for
    \code{method="profile"}? (See \code{signames} argument to
    \code{\link{profile}}).}
  \item{\dots}{additional parameters to be passed to
    \code{\link{profile.merMod}} or \code{\link{bootMer}}, respectively.}
}
\value{
  a numeric table of confidence intervals
}
\description{
  Compute confidence intervals on the parameters of a \code{*lmer()}
  model fit (of class\code{"\linkS4class{merMod}"}).
}
\details{
  Depending on the \code{method} specified, \code{confint()} computes
  confidence intervals by
  \describe{
    \item{\code{"profile"}:}{computing a
      likelihood profile and finding the appropriate cutoffs
      based on the likelihood ratio test;}
    \item{\code{"Wald"}:}{approximate
      the confidence intervals (of fixed-effect parameters
      only) based on the estimated local curvature of the
      likelihood surface;}
    \item{\code{"boot"}:}{perform parametric
      bootstrapping with confidence intervals computed from the
      bootstrap distribution according to \code{boot.type} (see
      \code{\link{boot.ci}}).}
  }
}
\examples{
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
fm1W <- confint(fm1, method="Wald")# very fast, but ....
fm1W
testLevel <- if (nzchar(s <- Sys.getenv("LME4_TEST_LEVEL"))) as.numeric(s) else 1
if(interactive() || testLevel >= 3) {
 ## ~20 seconds, MacBook Pro laptop
 system.time(fm1P <- confint(fm1, method="profile", ## default
                             oldNames = FALSE))
 ## ~ 40 seconds
 system.time(fm1B <- confint(fm1,method="boot",
                             .progress="txt", PBargs=list(style=3)))
} else 
  load(system.file("testdata","confint_ex.rda",package="lme4"))
fm1P
fm1B
}

