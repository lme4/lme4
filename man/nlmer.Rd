\name{nlmer}
\alias{nlmer}
\title{Fit Nonlinear Mixed-Effects Models}
\usage{
  nlmer(formula, data, control = list(), start = NULL,
    verbose = 0L, nAGQ = 1L, subset, weights, na.action,
    offset, contrasts = NULL, devFunOnly = 0L,
    tolPwrss = 1e-10,
    optimizer = c("NelderMead", "bobyqa"), ...)
}
\arguments{
  \item{formula}{a nonlinear mixed model formula (see
  detailed documentation)}

  \item{start}{starting estimates for the nonlinear model
  parameters, as a named numeric vector or as a list with
  components \describe{ \item{nlpars}{required numeric
  vector of starting values for the nonlinear model
  parameters} \item{theta}{optional numeric vector of
  starting values for the covariance parameters} }}

  \item{\dots}{other potential arguments.  A \code{method}
  argument was used in earlier versions of the package. Its
  functionality has been replaced by the \code{nAGQ}
  argument.}

  \item{data}{an optional data frame containing the
  variables named in \code{formula}.  By default the
  variables are taken from the environment from which
  \code{lmer} is called.}

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
  iterations in \code{glmer()}.  Default is
  \code{0.000001}.} }}

  \item{verbose}{integer scalar.  If \code{> 0} verbose
  output is generated during the optimization of the
  parameter estimates.  If \code{> 1} verbose output is
  generated during the individual PIRLS steps.}

  \item{nAGQ}{integer scalar - the number of points per
  axis for evaluating the adaptive Gauss-Hermite
  approximation to the log-likelihood.  Applies only to
  \code{glmer} and defaults to 1, corresponding to the
  Laplace approximation.  Values greater than 1 produce
  greater accuracy in the evaluation of the log-likelihood
  at the expense of speed.  A value of zero use a faster
  but less exact form of parameter estimation for GLMMs by
  optimizing the random effects and the fixed-effects
  coefficients in the penalized iteratively reweighted
  least squares step.}

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

  \item{tolPwrss}{numeric scalar - the tolerance for
  declaring convergence in the penalized iteratively
  weighted residual sum-of-squares step.  Defaults to
  1e-10.}

  \item{optimizer}{which optimizer to use in the second
  phase of optimization.  The default is
  \code{\link{NelderMead}} and the alternative is
  \code{\link{bobyqa}}.  For difficult model fits we have
  found \code{\link{NelderMead}} to be more reliable but
  occasionally slower than \code{\link{bobyqa}}.}
}
\description{
  Fit a nonlinear mixed-effects model
}
\details{
  Fit nonlinear mixed-effects models, such as those used in
  population pharmacokinetics.
}
\examples{
## nonlinear mixed models --- 3-part formulas ---

(nm1 <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree,
             Orange, start = c(Asym = 200, xmid = 725, scal = 350)))
(nm1a <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree,
              Orange, start = c(Asym = 200, xmid = 725, scal = 350),
              nAGQ = 0L))
}
\keyword{models}

