\name{profile-methods}
\title{Profile method for merMod objects}
\docType{methods}
\alias{as.data.frame.thpr}
\alias{log.thpr}
\alias{logProf}
\alias{profile-methods}
\alias{profile.merMod}
\description{
  Methods for \code{\link{profile}}() of [ng]\code{\link{lmer}} fitted
  models.

  The \code{log()} method and the more flexible \code{logProf()}
  utility transform a lmer profile in one where logarithms of standard deviations
  are used, see Details.
}
\usage{
\method{profile}{merMod}(fitted, which = 1:nptot, alphamax = 0.01,
	maxpts = 100, delta = cutoff/8, verbose = 0, devtol = 1e-09,
        maxmult = 10, startmethod = "prev", optimizer = "bobyqa",
	signames = TRUE, ...)
\method{as.data.frame}{thpr} (x, ...)
\method{log}{thpr}(x, base = exp(1))
logProf(x, base = exp(1), ranef = TRUE,
           sigIni = if(ranef) "sig" else "sigma")
}
\arguments{
  \item{fitted}{a fitted model, e.g., the result of \code{\link{lmer}(..)}.}
  \item{which}{integer or character vector indicating which parameters
    to profile: default is all parameters.  For integer, i.e., indexing,
    the parameters are ordered as follows:
    \describe{
      \item{(1)}{random effects (theta) parameters; these are ordered as
	in \code{getME(.,"theta")}, i.e., as the lower triangle of a
	matrix with standard deviations on the diagonal and correlations
	off the diagonal.}
      \item{(2)}{residual standard deviation (or scale parameter for GLMMs
        where appropriate);}
      \item{(3)}{fixed effect (beta) parameters.}
    }
    In addition, \code{which} may be a character, containing
    \code{"beta_"} or \code{"theta_"} which means the fixed and random
    effects parameters, respectively, or also containing parameter
    names, such as \code{".sigma"} or \code{"(Intercept)"}.
  }
  \item{alphamax}{a number in \eqn{(0,1)}, such that \code{1 - alphamax}
    is the maximum alpha value for likelihood ratio confidence
    regions; used to establish the range of values to be profiled.}
  \item{maxpts}{maximum number of points (in each direction, for each
    parameter) to evaluate in attempting to construct the profile.}
  \item{delta}{stepping scale for deciding on next point to profile.}
  \item{verbose}{level of output from internal calculations.}
  \item{devtol}{tolerance for fitted deviances less than
    baseline (supposedly minimum) deviance.}
  \item{maxmult}{maximum multiplier of the original step size allowed,
    defaults to 10.}
  \item{startmethod}{method for picking starting conditions for
    optimization (STUB).}
  \item{optimizer}{(character or function) optimizer to use (see
    \code{\link{lmer}} for details).}
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
}
\value{
  \code{profile(<merMod>)} returns
  an object of S3 class \code{"thpr"},% = th[eoretical]?? pr[ofile]
  \code{\link{data.frame}}-like.  Methods for such a profile object
  are notably \code{\link{confint}()} and the three plotting methods
  (which require the \pkg{lattice} package),
  \code{\link[=xyplot.thpr]{xyplot}}, \code{densityplot}, and
  \code{splom}.

  Further, \code{\link{log}()} (see above) and \code{\link{as.data.frame}()}.
}
\details{
  The \code{\link{log}} method and the more flexible \code{logProf()}
  function transform the profile into one where \eqn{\log(\sigma)} is
  used instead of \eqn{\sigma}.
  By default all sigmas, including the standard deviations of the random
  effects are transformed, i.e., it returns a profile with all the .sigNN
  parameters replaced by .lsigNN.  If \code{ranef} is false, only
  \code{".sigma"}, the standard deviation of the errors, is transformed
  as it should never be zero, whereas a random effect .sigNN is zero,
  when the corresponding random effect is entirely absent.
  \cr
  The forward and backward splines for the log-transformed parameters
  are recalculated.

  Methods for function \code{\link{profile}} (package
  \pkg{stats}), here for profiling (fitted) mixed effect models.

%% FIXME: ../inst/doc/profiling.txt  contains  motivation and more by
%% Doug Bates. Should add here (partly), or "link to there".
}
\seealso{
  The plotting methods \code{\link[=xyplot.thpr]{xyplot}} etc, for class
  \code{"thpr"}; \code{\link{varianceProf}} for transformation (from
  st.dev.) to variance scale.

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
stopifnot(all.equal(CIpr,
  array(c(12.1985292, 38.2299848, 1486.4515,
          84.0630513, 67.6576964, 1568.54849), dim = 3:2,
        dimnames = list(c(".sig01", ".sigma", "(Intercept)"),
                        c("2.5 \%", "97.5 \%"))),
                    tol= 1e-07))# 1.37e-9 {64b}
}
require(lattice)
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
 ## running time ~9 seconds on a modern machine:
 print( system.time(pr4 <- profile(gm1)) )
 print( xyplot(pr4,layout=c(5,1),as.table=TRUE) )
 print( xyplot(log(pr4), absVal=TRUE) ) # log(sigma_1)
 print( splom(pr4) )
 print( system.time( # quicker: only sig01 and one fixed effect
     pr2 <- profile(gm1, which=c("theta_", "period2"))))
 print( confint(pr2) )
}% doMore
}
\keyword{methods}
