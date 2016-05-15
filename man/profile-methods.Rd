\name{profile-methods}
\title{Profile method for merMod objects}
\docType{methods}
\alias{as.data.frame.thpr}
\alias{log.thpr}
\alias{logProf}
\alias{varianceProf}
\alias{profile-methods}
\alias{profile.merMod}
\description{
  Methods for \code{\link{profile}}() of [ng]\code{\link{lmer}} fitted
  models.

  The \code{log()} method and the more flexible \code{logProf()}
  utility transform a \code{lmer} profile into one where logarithms of standard deviations
  are used, while \code{varianceProf} converts from the
  standard-deviation to the variance scale; see Details.
}
\usage{
\method{profile}{merMod}(fitted, which = NULL, alphamax = 0.01,
	maxpts = 100, delta = NULL,
        delta.cutoff = 1/8, verbose = 0, devtol = 1e-09,
        maxmult = 10, startmethod = "prev", optimizer = NULL,
        control=NULL, signames = TRUE,
        parallel = c("no", "multicore", "snow"),
        ncpus = getOption("profile.ncpus", 1L), cl = NULL,
        prof.scale = c("sdcor","varcov"),
        \dots)
\method{as.data.frame}{thpr} (x, ...)
\method{log}{thpr}(x, base = exp(1))
logProf(x, base = exp(1), ranef = TRUE,
           sigIni = if(ranef) "sig" else "sigma")
varianceProf(x, ranef = TRUE)
}
\arguments{
  \item{fitted}{a fitted model, e.g., the result of \code{\link{lmer}(..)}.}
  \item{which}{NULL value,
    integer or character vector indicating which parameters
    to profile: default (NULL) is all parameters.  For integer, i.e., indexing,
    the parameters are ordered as follows:
    \describe{
      \item{(1)}{random effects (theta) parameters; these are ordered as
	in \code{getME(.,"theta")}, i.e., as the lower triangle of a
	matrix with standard deviations on the diagonal and correlations
	off the diagonal.}
      \item{(2)}{residual standard deviation (or scale parameter for GLMMs
        where appropriate).}
      \item{(3)}{fixed effect (beta) parameters.}
    }
    Alternatively, \code{which} may be a character, containing
    \code{"beta_"} or \code{"theta_"} denoting the fixed or random
    effects parameters, respectively, or also containing parameter
    names, such as \code{".sigma"} or \code{"(Intercept)"}.
  }
  \item{alphamax}{a number in \eqn{(0,1)}, such that \code{1 - alphamax}
    is the maximum alpha value for likelihood ratio confidence
    regions; used to establish the range of values to be profiled.}
  \item{maxpts}{maximum number of points (in each direction, for each
    parameter) to evaluate in attempting to construct the profile.}
  \item{delta}{stepping scale for deciding on next point to profile.
    The code uses the local derivative of the profile at the current
    step to establish a change in the focal parameter that will lead
    to a step of \code{delta} on the square-root-deviance scale.
    If \code{NULL}, the \code{delta.cutoff} parameter will be used
    to determine the stepping scale.}
  \item{delta.cutoff}{stepping scale (see \code{delta})
    expressed as a fraction of the
    target maximum value of the profile on the square-root-deviance
    scale.  Thus a \code{delta.cutoff} setting of \code{1/n} will
    lead to a profile with approximately \code{2*n} calculated points
    for each parameter (i.e., \code{n} points in each direction,
    below and above the estimate for each parameter).}
  \item{verbose}{level of output from internal calculations.}
  \item{devtol}{tolerance for fitted deviances less than
    baseline (supposedly minimum) deviance.}
  \item{maxmult}{maximum multiplier of the original step size allowed,
    defaults to 10.}
  \item{startmethod}{method for picking starting conditions for
    optimization (STUB).}
  \item{optimizer}{(character or function) optimizer to use (see
    \code{\link{lmer}} for details); default is to use the optimizer
    from the original model fit.}
  \item{control}{a \code{\link{list}} of options controlling the
    profiling (see \code{\link{lmerControl}}): default is to use the
    control settings from the original model fit.}
  \item{signames}{logical indicating if abbreviated names of the form
    \code{.sigNN} should be used; otherwise, names are more meaningful
    (but longer) of the form \code{(sd|cor)_(effects)|(group)}.  Note
    that some code for profile transformations (e.g., \code{log()} or
    \code{\link{varianceProf}}) depends on \code{signames==TRUE}.}
  \item{\dots}{potential further arguments for various methods.}
  \item{x}{an object of class \code{thpr} (i.e., output of
    \code{profile})}
  \item{base}{the base of the logarithm.  Defaults to natural
    logarithms.}
  \item{ranef}{logical indicating if the sigmas of the random effects
    should be \code{log()} transformed as well.  If false, only
    \eqn{\sigma} (standard deviation of errors) is transformed.}
  \item{sigIni}{character string specifying the initial part of the
    sigma parameters to be log transformed.}
  \item{parallel}{The type of parallel operation to be used (if any).
    If missing, the
    default is taken from the option \code{"profile.parallel"} (and if that
    is not set, \code{"no"}).}
  \item{ncpus}{integer: number of processes to be used in parallel operation:
    typically one would choose this to be the number of available CPUs.}
  \item{cl}{An optional \pkg{parallel} or \pkg{snow} cluster for use if
    \code{parallel = "snow"}.  If not supplied, a cluster on the
    local machine is created for the duration of the \code{profile}
    call.}
  \item{prof.scale}{whether to profile on the standard
    deviation-correlation scale (\code{"sdcor"}) or on
    the variance-covariance scale (\code{"varcov"})
  }
    
}
\value{
  \code{profile(<merMod>)} returns an object of S3 class
  \code{"thpr"}, %% = th[eta] pr[ofile], now a misnomer, as we also profile beta's
  which is \code{\link{data.frame}}-like.
  Notable methods for such a profile object
  \code{\link{confint}()}, which returns the
  confidence intervals based on the profile,
  and three plotting methods
  (which require the \pkg{lattice} package),
  \code{\link[=xyplot.thpr]{xyplot}}, \code{densityplot}, and
  \code{splom}.

  In addition, the
  \code{\link{log}()} (see above) and \code{\link{as.data.frame}()}
  methods can transform \code{"thpr"} objects in useful ways.
}
\details{
  The \code{\link{log}} method and the more flexible \code{logProf()}
  function transform the profile into one where \eqn{\log(\sigma)} is
  used instead of \eqn{\sigma}.
  By default all sigmas including the standard deviations of the random
  effects are transformed i.e., the methods return a profile with all
  of the \code{.sigNN}
  parameters replaced by \code{.lsigNN}.  If \code{ranef} is false, only
  \code{".sigma"}, the standard deviation of the errors, is transformed
  (as it should never be zero, whereas random effect standard
  deviations (\code{.sigNN}) can be reasonably be zero).
  \cr
  The forward and backward splines for the log-transformed parameters
  are recalculated.
  Note that correlation parameters are not handled sensibly at present
  (i.e., they are logged rather than taking a more applicable
  transformation such as an arc-hyperbolic tangent,
  \code{atanh(x)}=\eqn{\log((1+x)/(1-x))/2}{log((1+x)/(1-x))/2}).
  
  The \code{varianceProf} function works similarly, including
  non-sensibility for correlation parameters, by squaring all
  parameter values, changing the names by appending \code{sq}
  appropriately (e.g. \code{.sigNN} to \code{.sigsqNN}).
  Setting \code{prof.scale="varcov"} in the original
  \code{profile()} call is a more computationally
  intensive, but more correct, way to compute confidence
  intervals for covariance parameters.

  Methods for function \code{\link{profile}} (package
  \pkg{stats}), here for profiling (fitted) mixed effect models.

%% FIXME: ../inst/doc/profiling.txt  contains  motivation and more by
%% Doug Bates. Should add here (partly), or "link to there".
}
\seealso{
  The plotting methods \code{\link[=xyplot.thpr]{xyplot}} etc, for class
  \code{"thpr"}.

  For (more expensive) alternative confidence intervals:
  \code{\link{bootMer}}.
}
\examples{
fm01ML <- lmer(Yield ~ 1|Batch, Dyestuff, REML = FALSE)
system.time(
  tpr  <- profile(fm01ML, optimizer="Nelder_Mead", which="beta_")
)## fast; as only *one* beta parameter is profiled over
## full profiling (default which means 'all) needs
## ~2.6s (on a 2010 Macbook Pro)
system.time( tpr  <- profile(fm01ML))
## ~1s, + possible warning about bobyqa convergence
(confint(tpr) -> CIpr)
\donttest{% too much precision (etc). but just FYI:
stopifnot(all.equal(unname(CIpr),
  array(c(12.1985292, 38.2299848, 1486.4515,
          84.0630513, 67.6576964, 1568.54849), dim = 3:2),
                    tol= 1e-07))# 1.37e-9 {64b}
}
library("lattice")
xyplot(tpr)
xyplot(tpr, absVal=TRUE) # easier to see conf.int.s (and check symmetry)
xyplot(tpr, conf = c(0.95, 0.99), # (instead of all five 50, 80,...)
       main = "95\% and 99\% profile() intervals")
xyplot(logProf(tpr, ranef=FALSE),
       main = expression("lmer profile()s"~~ log(sigma)*" (only log)"))
densityplot(tpr, main="densityplot( profile(lmer(..)) )")
densityplot(varianceProf(tpr), main=" varianceProf( profile(lmer(..)) )")
splom(tpr)
splom(logProf(tpr, ranef=FALSE))
doMore <- lme4:::testLevel() > 1 %% even more --> ../tests/profile.R
if(doMore) { ## not typically, for time constraint reasons
 ## Batch and residual variance only
 system.time(tpr2 <- profile(fm01ML, which=1:2, optimizer="Nelder_Mead"))
 print( xyplot(tpr2) )
 print( xyplot(log(tpr2)) )# log(sigma) is better
 print( xyplot(logProf(tpr2, ranef=FALSE)) )

 ## GLMM example
 gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
	     data = cbpp, family = binomial)
 ## running ~ 10-12 seconds on a modern machine {-> "verbose" while you wait}:
 print( system.time(pr4 <- profile(gm1, verbose=TRUE)) )
 print( xyplot(pr4, layout=c(5,1), as.table=TRUE) )
 print( xyplot(log(pr4), absVal=TRUE) ) # log(sigma_1)
 print( splom(pr4) )
 print( system.time( # quicker: only sig01 and one fixed effect
     pr2 <- profile(gm1, which=c("theta_", "period2"))))
 print( confint(pr2) )
 ## delta..: higher underlying resolution, only for 'sigma_1':
 print( system.time(
     pr4.hr <- profile(gm1, which="theta_", delta.cutoff=1/16)))
 print( xyplot(pr4.hr) )
}% doMore
}
\keyword{methods}
