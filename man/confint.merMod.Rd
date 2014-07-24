\name{confint.merMod}
\alias{confint.merMod}
\title{Compute Confidence Intervals for Parameters of a [ng]lmer Fit}
\description{
  Compute confidence intervals on the parameters of a \code{*lmer()}
  model fit (of class\code{"\linkS4class{merMod}"}).
}
\usage{
\method{confint}{merMod}(object, parm, level = 0.95,
	method = c("profile", "Wald", "boot"), zeta,
	nsim = 500,
        boot.type = c("perc","basic","norm"), quiet = FALSE,
	oldNames = TRUE, ...)
}
\arguments{
  \item{object}{a fitted [ng]lmer model}
  \item{parm}{parameters of which intervals are sought. Specified by
    integer position, or for \code{method = "profile"}, e.g., via a
    \code{\link{character}} vector, see Note below, and
    \code{\link[=profile.merMod]{profile}}.}
  \item{level}{confidence level \eqn{< 1}, typically above 0.90.}
  \item{method}{a \code{\link{character}} string determining the method
    for computing the confidence intervals.}
  \item{zeta}{(for \code{method = "profile"} only:) likelihood cutoff
    (if not specified, as by default, computed from \code{level}).}
  \item{nsim}{number of simulations for parametric bootstrap intervals.}
  \item{boot.type}{bootstrap confidence interval type, as described
    in \code{\link{boot.ci}}. (Methods \sQuote{stud} and \sQuote{bca}
    are unavailable because they require additional components to be
    calculated.)}
  \item{quiet}{(logical) suppress messages about computationally intensive profiling?}
  \item{oldNames}{(logical) use old-style names for
    \code{method="profile"}? (See \code{signames} argument to
    \code{\link{profile}}).}
  \item{\dots}{additional parameters to be passed to
    \code{\link{profile.merMod}} or \code{\link{bootMer}}, respectively.}
}
\value{
  a numeric table (\code{\link{matrix}} with column and row names) of
  confidence intervals.
}
\note{
  The default method \code{"profile"} amounts to \preformatted{
    confint(profile(object, which=parm), signames=oldNames, ...),
            level, zeta)}
  where the \code{\link{profile}} method \code{\link{profile.merMod}}
  does almost all the computations.  Therefore it is typically
  advisable to store the \link[=profile.merMod]{profile(.)} result, say
  in \code{pp}, and then use \code{confint(pp, level=*)} e.g., for
  different levels.
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
      \code{\link[boot]{boot.ci}}).}
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

