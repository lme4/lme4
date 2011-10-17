\name{lmer}
\alias{lmer}
\alias{glmer}
\alias{nlmer}
\concept{LMM}% << so it's found; ditto for these:
\concept{GLMM}
\concept{NLMM}
\title{Fit Mixed-Effects Models}
\description{
  Fit a linear mixed model or a generalized linear mixed model or a
  nonlinear mixed model.
}
\usage{
lmer(formula, data, REML = TRUE, sparseX = FALSE,
     control = list(), start = NULL, verbose = 0L,
     doFit = TRUE, %compDev = TRUE, _support too - FIXME? _
     subset, weights,
     na.action, offset, contrasts = NULL,
     devFunOnly=FALSE,
     \dots)

glmer(formula, data, family = gaussian, sparseX = FALSE,
      control = list(), start = NULL, verbose = 0, nAGQ = 1,
      doFit = TRUE, compDev = TRUE, subset, weights, na.action, offset,
      contrasts = NULL, mustart, etastart, devFunOnly = FALSE,
      tolPwrss = 1e-06, \dots)

nlmer(formula, data, family = gaussian, start = NULL,
      verbose = 0L, nAGQ = 1L, doFit = TRUE, subset, weights,
      na.action, offset, contrasts = NULL, devFunOnly = FALSE, \dots)
}
\arguments{
  \item{formula}{a two-sided linear formula object describing the
    fixed-effects part of the model, with the response on the left of a
    \code{~} operator and the terms, separated by \code{+} operators, on
    the right.  The vertical bar character \code{"|"} separates an
    expression for a model matrix and a grouping factor.}
  \item{data}{an optional data frame containing the variables named in
    \code{formula}.  By default the variables are taken from the
    environment from which \code{lmer} is called.}
  \item{sparseX}{logical - should a sparse model matrix be used for the
    fixed-effects terms?  Defaults to \code{FALSE}.}
  \item{family}{a GLM family, see \code{\link[stats]{glm}} and
    \code{\link[stats]{family}}.}
  \item{REML}{logical - Should the estimates
    be chosen to optimize the REML criterion (as opposed to the
    log-likelihood)?  Defaults to \code{TRUE}.}
  \item{nAGQ}{an integer - the number of points per axis for
    evaluating the adaptive Gauss-Hermite approximation to the
    log-likelihood.  This defaults to 1, corresponding to the Laplacian
    approximation.  Values greater than 1 produce greater accuracy in
    the evaluation of the log-likelihood at the expense of speed.  A
    value of zero use a quicker but less exact form of parameter
    estimation for nonlinear models by optimizing the random effects and
    the fixed-effects coefficients in the penalized least squares
    step.}
  \item{control}{a list of control parameters. See below for details.}
  \item{start}{a named list of starting values for the parameters in the
    model.  If the list
    contains components named \code{fixef} and/or \code{theta}, these are
    used as the starting values for those slots.  (Setting starting
    values for \code{fixef} has no effect for a linear mixed model
    because the fixed-effects parameters do not appear in the profiled
    deviance.)  In \code{lmer} and \code{glmer} a numeric
    \code{start} argument of the appropriate length is used as the
    starting value of \code{theta}.  In \code{nlmer} a numeric
      \code{start} argument is used as the starting values of the
      \code{fixef} slot.}
  \item{doFit}{logical scalar. When \code{doFit = FALSE} the model is
    not fit but instead a structure with the model matrices for the
    random-effects terms is returned, so they can be modified for
    special model forms. When \code{doFit = TRUE}, the default, the
    model is fit immediately.}
  \item{compDev}{logical - should compiled code be used for the deviance
    evaluation during the optimization of the parameters estimates?}
  \item{subset, weights, na.action, offset, contrasts}{further model
    specification arguments as in \code{\link[stats]{lm}}; see there for
    details.}
  \item{mustart, etastart}{optional starting values on the scale of the
    conditional mean or the unbounded predictor as in
    \code{\link[stats]{glm}}; see there for details.}
  \item{verbose}{integer scalar.  If \code{> 0} verbose output is
    generated during the optimization of the parameter estimates.  If
    \code{> 1} verbose output is generated during the individual PIRLS
    steps.}
  \item{devFunOnly}{logical - return only the deviance evaluation
    function - useful mainly for debugging.}% but MM thinks we want to keep.
  \item{tolPwrss}{numeric scalar - the tolerance for declaring
      convergence in the penalized weighted residual sum-of-squares
      step.  Defaults to 1e-06.}
  \item{\dots}{other potential arguments.  A \code{method} argument was
    used in earlier versions of the package. Its functionality has been
    replaced by the \code{REML} and \code{nAGQ} arguments.}
}
\details{
  The \code{lmer} and \code{glmer} functions are nearly interchangeable.
  If \code{lmer} is called with a non-default \code{family} argument the
  call is replaced by a call to \code{glmer} with the current arguments.
  If \code{glmer} is called with the default \code{family}, namely the
  \code{\link[stats:family]{gaussian}} family with the identity link, then the call is
  replaced by a call to \code{lmer} with the current arguments.  (They
  are described as \dQuote{nearly} interchangeable because the
  \code{REML} argument only applies to calls to \code{lmer} and the
  \code{nAGQ} argument only applies to calls to \code{glmer}.)

  % The \code{lmer2} function has the same purpose as \code{lmer} but uses
  % objects from reference classes defined in an Rcpp Module.

  Additional standard arguments to model-fitting functions can be passed
  to \code{lmer} or \code{glmer} or \code{nlmer}:
  \describe{
    \item{subset}{an optional expression indicating the subset of the rows of
      \code{data} that should be used in the fit. This can be a logical
      vector, or a numeric vector indicating which observation numbers are
      to be included, or a  character  vector of the row names to be
      included.  All observations are included by default.}
    \item{na.action}{a function that indicates what should happen when the
      data contain \code{NA}s.  The default action (\code{na.fail}) prints
      an error message and terminate if there are any incomplete
      observations.}
    \item{control}{a named list of control parameters for the estimation
      algorithm, specifying only the ones to be changed from their
      default values.  Hence defaults to an empty list.\cr
      Possible control options and their default values are:
      \describe{
	\item{\code{msVerbose}:}{a logical value passed as the
	  \code{trace} argument to \code{nlminb} (see documentation on
	  that function).  Default is \code{getOption("verbose")}.}
	\item{\code{maxIter}:}{a positive integer passed as the
	  \code{maxIter} argument to \code{nlminb} (see documentation on
	  that function).  Default is \code{300}.}
	\item{\code{maxFN}:}{a positive integer specifying the
	  maximum number of evaluations of the deviance function allowed
	  during the optimization. Default is \code{900}}.
	\item{\code{tol}:}{a positive number specifying the
	  convergence tolerance, currently only for the PWRSS iterations
	  in \code{glmer()}.  Default is \code{0.000001}}.
      }
    }
    \item{model, x}{logicals.  If \code{TRUE} the corresponding
      components of the fit (the model frame, the model matrices)
      are returned.}
  }
}
\value{
  An object of class \code{"\linkS4class{merMod}"}, for which many methods
  are available.  See there for details.
}
\seealso{
  The \code{\linkS4class{merMod}} class, \code{\link[stats]{lm}}
}
\examples{
## linear mixed models - reference values from older code
(fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
(fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))
anova(fm1, fm2)

## dense vs sparse X  --------------------------
## use more sensible example !!
fm3.d <- lmer(Yield ~ 1|Batch, Dyestuff2)
fm3.s <- lmer(Yield ~ 1|Batch, Dyestuff2, sparseX=TRUE)#-> warning
## check "equality"
stopifnot(all.equal( coef(fm3.d), coef(fm3.s), tol = 1e-14),
          all.equal(sigma(fm3.d),sigma(fm3.s), tol = 1e-14),
          TRUE)


## generalized linear mixed model
(gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
              data = cbpp, family = binomial))
## using nAGQ=0L only gets close to the optimum
(gm1a <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
               cbpp, binomial, nAGQ = 0L))

## GLMM with individual-level variability (accounting for overdispersion)
cbpp$obs <- 1:nrow(cbpp)
(gm2 <- glmer(cbind(incidence, size - incidence) ~ period +
    (1 | herd) +  (1|obs),
              family = binomial, data = cbpp))
anova(gm1,gm2)


## nonlinear mixed models --- 3-part formulas ---

nm1 <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~
             0 + Asym + xmid + scal + (0 + Asym|Tree),
             Orange, start = c(Asym = 200, xmid = 725, scal = 350),
             verbose = 1)
nm1
nm1a <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~
              0 + Asym + xmid + scal + (0 + Asym|Tree),
              Orange, start = c(Asym = 200, xmid = 725, scal = 350),
              verbose = 1, nAGQ = 0L)
nm1a
}
\keyword{models}
\keyword{methods}
