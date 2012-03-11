\name{lmer}
\alias{lmer}
\title{Fit Linear Mixed-Effects Models}
\usage{
  lmer(formula, data, REML = TRUE, sparseX = FALSE,
    control = list(), start = NULL, verbose = 0L, subset,
    weights, na.action, offset, contrasts = NULL,
    devFunOnly = FALSE, optimizer="Nelder_Mead", ...)
}
\arguments{
  \item{formula}{a two-sided linear formula object
  describing the fixed-effects part of the model, with the
  response on the left of a \code{~} operator and the
  terms, separated by \code{+} operators, on the right.
  The vertical bar character \code{"|"} separates an
  expression for a model matrix and a grouping factor.}

  \item{data}{an optional data frame containing the
  variables named in \code{formula}.  By default the
  variables are taken from the environment from which
  \code{lmer} is called.}

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
    evaluation function.}
  \item{optimizer}{optimizer to use}

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
\examples{
## linear mixed models - reference values from older code
(fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
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

