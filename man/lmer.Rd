\name{lmer}
\alias{lmer}
\title{Fit Linear Mixed-Effects Models}
\usage{
  lmer(formula, data = NULL, REML = TRUE, sparseX = FALSE,
    control = list(), start = NULL, verbose = 0L, subset,
    weights, na.action, offset, contrasts = NULL,
    devFunOnly = FALSE, optimizer = "Nelder_Mead", ...)
}
\arguments{
  \item{formula}{a two-sided linear formula object
  describing both the fixed-effects and fixed-effects part
  of the model, with the response on the left of a \code{~}
  operator and the terms, separated by \code{+} operators,
  on the right.  Random-effects terms are distinguished by
  vertical bars (\code{"|"}) separating expressions for
  design matrices from grouping factors.}

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

  \item{REML}{logical scalar - Should the estimates be
  chosen to optimize the REML criterion (as opposed to the
  log-likelihood)?  Defaults to \code{TRUE}.}

  \item{sparseX}{logical - should a sparse model matrix be
  used for the fixed-effects terms?  Defaults to
  \code{FALSE}. Currently inactive.}

  \item{control}{a named list of control parameters for the
  estimation algorithm, specifying only the ones to be
  changed from their default values.  Hence defaults to an
  empty list.\cr Possible control options and their default
  values are: \describe{ \item{\code{msVerbose}:}{a logical
  value passed as the \code{trace} argument to
  \code{nlminb} (see documentation on that function).
  Default is \code{getOption("verbose")}.}
  \item{\code{maxIter}:}{a positive integer passed as the
  \code{maxIter} argument to \code{nlminb} (see
  documentation on that function).  Default is \code{300}.}
  \item{\code{maxFN}:}{a positive integer specifying the
  maximum number of evaluations of the deviance function
  allowed during the optimization. Default is \code{900}.}
  \item{\code{tol}:}{a positive number specifying the
  convergence tolerance, currently only for the PWRSS
  iterations in \code{\link{glmer}}.  Default is
  \code{0.000001}.} }}

  \item{start}{a named list of starting values for the
  parameters in the model.  For \code{lmer} this can be a
  numeric vector or a list with one component named
  \code{"theta"}. Infrequently used.}

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

  \item{weights}{an optional vector of \sQuote{prior
  weights} to be used in the fitting process.  Should be
  \code{NULL} or a numeric vector.}

  \item{na.action}{a function that indicates what should
  happen when the data contain \code{NA}s.  The default
  action (\code{na.fail}) prints an error message and
  terminates if there are any incomplete observations.}

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

  \item{optimizer}{character - name of optimizing function.
  The built-in optimizers are \code{\link{Nelder_Mead}} and
  \code{\link[minqa]{bobyqa}} (from the \pkg{minqa}
  package. Any minimizing function that allows box
  constraints can be used provided that it (1) takes input
  parameters \code{fn} (function to be optimized),
  \code{par} (starting parameter values), \code{lower}
  (lower bounds) and \code{control} (control parameters,
  passed through from the \code{control} argument) and (2)
  returns a list with (at least) elements \code{par}
  (best-fit parameters), \code{fval} (best-fit function
  value), \code{conv} (convergence code) and (optionally)
  \code{message} (informational message, or explanation of
  convergence failure).  Special provisions are made for
  \code{\link{bobyqa}}, \code{\link{Nelder_Mead}}, and
  optimizers wrapped in the \pkg{optimx} package; to use
  \pkg{optimx} optimizers (including \code{L-BFGS-B} from
  base \code{\link{optim}} and \code{\link{nlminb}}), pass
  the \code{method} argument to \code{optim} in the
  \code{control} argument.}

  \item{\dots}{other potential arguments.  A \code{method}
  argument was used in earlier versions of the package. Its
  functionality has been replaced by the \code{REML}
  argument.}
}
\value{
  An object of class \code{"\linkS4class{merMod}"}, for
  which many methods are available.  See there for details.
}
\description{
  Fit a linear mixed model (LMM)
}
\details{
  \itemize{ \item{If the \code{formula} argument is
  specified as a character vector, the function will
  attempt to coerce it to a formula. However, this is not
  recommended (users who want to construct formulas by
  pasting together components are advised to use
  \code{\link{as.formula}}); model fits will work but
  subsequent methods such as \code{\link{drop1}},
  \code{\link{update}} may fail.} \item{Unlike some simpler
  modeling frameworks such as \code{\link{lm}} and
  \code{\link{glm}} which automatically detect perfectly
  collinear predictor variables, \code{[gn]lmer} cannot
  handle design matrices of less than full rank.  For
  example, in cases of models with interactions that have
  unobserved combinations of levels, it is up to the user
  to define a new variable (for example creating \code{ab}
  within the data from the results of
  \code{droplevels(interaction(a,b))}). } \item{the
  deviance function returned when \code{devFunOnly} is
  \code{TRUE} takes a single numeric vector argument,
  representing the \code{theta} vector. This vector defines
  the scaled variance-covariance matrices of the random
  effects, in the Cholesky parameterization.  For models
  with only simple (intercept-only) random effects,
  \code{theta} is a vector of the standard deviations of
  the random effects. For more complex or multiple random
  effects, running \code{getME(.,"theta")} to retrieve the
  \code{theta} vector for a fitted model and examining the
  names of the vector is probably the easiest way to
  determine the correspondence between the elements of the
  \code{theta} vector and elements of the lower triangles
  of the Cholesky factors of the random effects.} }
}
\examples{
## linear mixed models - reference values from older code
(fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
fm1_ML <- update(fm1,REML=FALSE)
(fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))
anova(fm1, fm2)
}
\seealso{
  The \code{\linkS4class{merMod}} class,
  \code{\link[stats]{lm}}
}
\concept{
  LMM
}
\keyword{models}

