\name{lmer}
\title{Fit Linear Mixed-Effects Models}
\alias{lmer}
\concept{ LMM }
\description{
  Fit a linear mixed-effects model (LMM) to data.
}
\usage{
  lmer(formula, data = NULL, REML = TRUE,
    control = lmerControl(), start = NULL, verbose = 0L,
    subset, weights, na.action, offset, contrasts = NULL,
    devFunOnly = FALSE, ...)
}
\arguments{
  \item{formula}{a two-sided linear formula object
  describing both the fixed-effects and random-effects part
  of the model, with the response on the left of a \code{~}
  operator and the terms, separated by \code{+} operators,
  on the right.  Random-effects terms are distinguished by
  vertical bars (\code{"|"}) separating expressions for
  design matrices from grouping factors. Two vertical bars
  (\code{"||"}) can be used to specify  multiple uncorrelated
  random effects for the same grouping variable.}


  \item{data}{an optional data frame containing the
  variables named in \code{formula}.  By default the
  variables are taken from the environment from which
  \code{lmer} is called. While \code{data} is optional, the
  package authors \emph{strongly} recommend its use,
  especially when later applying methods such as
  \code{update} and \code{drop1} to the fitted model
  (\emph{such methods are not guaranteed to work properly
  if \code{data} is omitted}). If \code{data} is omitted,
  variables will be taken from the environment of
  \code{formula} (if specified as a formula) or from the
  parent frame (if specified as a character vector).}

  \item{REML}{logical scalar - Should the estimates be chosen to
    optimize the REML criterion (as opposed to the log-likelihood)?}

  \item{control}{a list (of correct class, resulting from
    \code{\link{lmerControl}()} or \code{\link{glmerControl}()}
    respectively) containing control parameters, including the nonlinear
    optimizer to be used and parameters to be passed through to the
    nonlinear optimizer, see the \code{*lmerControl} documentation for
    details.}

  \item{start}{a named \code{\link{list}} of starting values for the
    parameters in the model.  For \code{lmer} this can be a numeric
    vector or a list with one component named \code{"theta"}.}

  \item{verbose}{integer scalar.  If \code{> 0} verbose
  output is generated during the optimization of the
  parameter estimates.  If \code{> 1} verbose output is
  generated during the individual PIRLS steps.}

  \item{subset}{an optional expression indicating the
  subset of the rows of \code{data} that should be used in
  the fit. This can be a logical vector, or a numeric
  vector indicating which observation numbers are to be
  included, or a character vector of the row names to be
  included.  All observations are included by default.}

  \item{weights}{an optional vector of \sQuote{prior weights} to be used
  in the fitting process.  Should be \code{NULL} or a numeric vector.
  Prior \code{weights} are \emph{not} normalized or standardized in any
  way.  In particular, the diagonal of the residual covariance matrix is
  the squared residual standard deviation parameter \code{\link{sigma}}
  times the vector of inverse \code{weights}.  Therefore, if the
  \code{weights} have relatively large magnitudes, then in order to
  compensate, the \code{\link{sigma}} parameter will also need to have a
  relatively large magnitude.}

  \item{na.action}{a function that indicates what should
  happen when the data contain \code{NA}s.  The default
  action (\code{na.omit}, inherited from the 'factory
  fresh' value of \code{getOption("na.action")}) strips any
  observations with any missing values in any variables.}

  \item{offset}{this can be used to specify an \emph{a
  priori} known component to be included in the linear
  predictor during fitting. This should be \code{NULL} or a
  numeric vector of length equal to the number of cases.
  One or more \code{\link{offset}} terms can be included in
  the formula instead or as well, and if more than one is
  specified their sum is used.  See
  \code{\link{model.offset}}.}

  \item{contrasts}{an optional list. See the
  \code{contrasts.arg} of \code{model.matrix.default}.}

  \item{devFunOnly}{logical - return only the deviance
  evaluation function. Note that because the deviance
  function operates on variables stored in its environment,
  it may not return \emph{exactly} the same values on
  subsequent calls (but the results should always be within
  machine tolerance).}

  \item{\dots}{other potential arguments.  A \code{method}
  argument was used in earlier versions of the package. Its
  functionality has been replaced by the \code{REML}
  argument.}
}
\value{
  An object of class \code{merMod}, for which many methods
  are available (e.g. \code{methods(class="merMod")})
}
\details{
  \itemize{ \item{If the \code{formula} argument is
  specified as a character vector, the function will
  attempt to coerce it to a formula. However, this is not
  recommended (users who want to construct formulas by
  pasting together components are advised to use
  \code{\link{as.formula}} or \code{\link{reformulate}});
  model fits will work but subsequent methods such as
  \code{\link{drop1}}, \code{\link{update}} may fail.}
  \item{Unlike some simpler modeling frameworks such as
  \code{\link{lm}} and \code{\link{glm}} which
  automatically detect perfectly collinear predictor
  variables, \code{[gn]lmer} cannot handle design matrices
  of less than full rank.  For example, in cases of models
  with interactions that have unobserved combinations of
  levels, it is up to the user to define a new variable
  (for example creating \code{ab} within the data from the
  results of \code{interaction(a,b,drop=TRUE)}). }
  \item{the deviance function returned when
  \code{devFunOnly} is \code{TRUE} takes a single numeric
  vector argument, representing the \code{theta} vector.
  This vector defines the scaled variance-covariance
  matrices of the random effects, in the Cholesky
  parameterization.  For models with only simple
  (intercept-only) random effects, \code{theta} is a vector
  of the standard deviations of the random effects. For
  more complex or multiple random effects, running
  \code{getME(.,"theta")} to retrieve the \code{theta}
  vector for a fitted model and examining the names of the
  vector is probably the easiest way to determine the
  correspondence between the elements of the \code{theta}
  vector and elements of the lower triangles of the
  Cholesky factors of the random effects.} }
}
\seealso{
  \code{\link[stats]{lm}}
}
\examples{
## linear mixed models - reference values from older code
(fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy))
summary(fm1)# (with its own print method)

str(terms(fm1))
stopifnot(identical(terms(fm1, fixed.only=FALSE),
                    terms(model.frame(fm1))))
attr(terms(fm1, FALSE), "dataClasses") # fixed.only=FALSE needed for dataCl.

fm1_ML <- update(fm1,REML=FALSE)
(fm2 <- lmer(Reaction ~ Days + (Days || Subject), sleepstudy))
anova(fm1, fm2)
sm2 <- summary(fm2)
print(fm2, digits=7, ranef.comp="Var") # the print.merMod()         method
print(sm2, digits=3, corr=FALSE)       # the print.summary.merMod() method

(vv <- vcov.merMod(fm2, corr=TRUE))
as(vv, "corMatrix")# extracts the ("hidden") 'correlation' entry in @factors

## Fit sex-specific variances by constructing numeric dummy variables
## for sex and sex:age; in this case the estimated variance differences
## between groups in both intercept and slope are zero ...
data(Orthodont,package="nlme")
Orthodont$nsex <- as.numeric(Orthodont$Sex=="Male")
Orthodont$nsexage <- with(Orthodont, nsex*age)
lmer(distance ~ age + (age|Subject) + (0+nsex|Subject) +
     (0 + nsexage|Subject), data=Orthodont)

}
\keyword{models}
