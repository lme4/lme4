\name{lmer}
\alias{lmer}
\alias{lmer2}
\alias{glmer}
\alias{nlmer}
\concept{GLMM}% << so it's found
\concept{NLMM}% << so it's found
\title{Fit Mixed-Effects Models}
\description{
  Fit a linear mixed model or a generalized linear mixed model or a
  nonlinear mixed model.
}
\usage{
lmer(formula, data, family = NULL, REML = TRUE,
     control = list(), start = NULL, verbose = FALSE,
     doFit = TRUE, subset, weights, na.action, offset,
     contrasts = NULL, model = TRUE, x = TRUE, \dots)
lmer2(formula, data, family = NULL, REML = TRUE,
      control = list(), start = NULL, verbose = FALSE,
      subset, weights, na.action, offset,
      contrasts = NULL, model = TRUE, x = TRUE, \dots)
glmer(formula, data, family = gaussian, start = NULL,
      verbose = FALSE, nAGQ = 1, doFit = TRUE, subset, weights,
      na.action, offset, contrasts = NULL, model = TRUE,
      control = list(), \dots)
nlmer(formula, data, start = NULL, verbose = FALSE, nAGQ = 1,
      doFit = TRUE, subset, weights, na.action,
      contrasts = NULL, model = TRUE, control = list(), \dots)
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
  \item{family}{
    a GLM family, see \code{\link[stats]{glm}} and
    \code{\link[stats]{family}}. If \code{family} is missing then a
    linear mixed model is fit; otherwise a generalized linear mixed
    model is fit.}
  \item{REML}{logical argument to \code{lmer} only. Should the estimates
    be chosen to optimize the REML criterion (as opposed to the
    log-likelihood)?  Defaults to \code{TRUE}.}
  \item{nAGQ}{a positive integer - the number of points per axis for
    evaluating the adaptive Gauss-Hermite approximation to the
    log-likelihood.  This defaults to 1, corresponding to the Laplacian
    approximation.  Values greater than 1 produce greater accuracy in
    the evaluation of the log-likelihood at the expense of speed.}
  \item{control}{a list of control parameters. See below for details.}
  \item{start}{a named list of starting values for the parameters in the
    model. If the list is of the same form as the \code{ST} slot, it is
    becomes the starting values of the \code{ST} slot.  If the list
    contains components named \code{fixef} and/or \code{ST}, these are
    used as the starting values for those slots.  (Setting starting
    values for \code{fixef} has no effect for a linear mixed model
    because the fixed-effects parameters do not appear in the profiled
    deviance.)  In \code{lmer} and \code{glmer} a numeric
    \code{start} argument of the appropriate length is used as the
    starting value of the parameter vector that determines the \code{ST}
    slot. In \code{nlmer} a numeric \code{start} argument is used as the
    starting values of the \code{fixef} slot.}
  \item{doFit}{logical scalar. When \code{doFit = FALSE} the model is
    not fit but instead a structure with the model matrices for the
    random-effects terms is returned, so they can be modified for
    special model forms. When \code{doFit = TRUE}, the default, the
    model is fit immediately.}
  \item{subset, weights, na.action, offset, contrasts}{further model
    specification arguments as in \code{\link[stats]{lm}}; see there for
    details.}
  \item{model}{logical scalar.  If \code{FALSE} the model frame in
    slot \code{frame} is truncated to zero rows.}
  \item{x}{logical scalar.  If \code{FALSE} the model matrix in
    slot \code{X} is truncated to zero rows.}
  \item{verbose}{logical scalar.  If \code{TRUE} verbose output is
    generated during the optimization of the parameter estimates.}
  \item{\dots}{other potential arguments.  A \code{method} argument was
    used in earlier versions of the package. Its functionality has been
    replaced by the \code{REML} and \code{nAGQ} arguments.}
}
\details{

  The \code{lmer} and \code{glmer} functions are nearly interchangeable.
  If \code{lmer} is called with a non-default \code{family} argument the
  call is replaced by a call to \code{glmer} with the current arguments.
  If \code{glmer} is called with the default \code{family}, namely the
  \code{\link[stats:family]{gaussian}} family with the identity link,
  then the call is replaced by a call to \code{lmer} with the current
  arguments.  (They are described as \dQuote{nearly} interchangeable
  because the \code{REML} argument only applies to calls to \code{lmer}
  and the \code{nAGQ} argument only applies to calls to \code{glmer}.)

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
      }
    }
    \item{model, x}{logicals.  If \code{TRUE} the corresponding
      components of the fit (the model frame, the model matrices)
      are returned.}
  }

  The \code{lmer2} name exists only for backwards compatibility.
  Calling this function simply produces an equivalent call to \code{lmer}.
}
\value{
  An object of class \code{"\linkS4class{mer}"}, for which many methods
  are available.  See there for details.
}
\seealso{
  The \code{\linkS4class{mer}} class, \code{\link[stats]{lm}}
}
\examples{
## linear mixed models
(fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
(fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))
anova(fm1, fm2)
## generalized linear mixed model
(gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
              family = binomial, data = cbpp))
## GLMM with individual-level variability (accounting for overdispersion)
cbpp$obs <- 1:nrow(cbpp)
(gm2 <- glmer(cbind(incidence, size - incidence) ~ period +
    (1 | herd) +  (1|obs),
              family = binomial, data = cbpp))
anova(gm1,gm2)
## nonlinear mixed models
(nm1 <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree,
              Orange, start = c(Asym = 200, xmid = 725, scal = 350)))
}
\keyword{models}
\keyword{methods}
