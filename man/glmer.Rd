\name{glmer}
\alias{glmer}
\title{Fit Generalized Linear Mixed-Effects Models}
\usage{
  glmer(formula, data, family = gaussian, sparseX = FALSE,
    control = list(), start = NULL, verbose = 0L,
    nAGQ = 1L, compDev = TRUE, subset, weights, na.action,
    offset, contrasts = NULL, mustart, etastart,
    devFunOnly = FALSE, tolPwrss = 1e-10,
    optimizer = c("NelderMead", "bobyqa"), ...)
}
\arguments{
  \item{family}{a GLM family, see \code{\link[stats]{glm}}
  and \code{\link[stats]{family}}.}

  \item{compDev}{logical scalar - should compiled code be
  used for the deviance evaluation during the optimization
  of the parameter estimates?  Defaults to \code{TRUE}.}

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

  \item{start}{a named list of starting values for the
  parameters in the model.  If the list contains components
  named \code{fixef} and/or \code{theta}, these are used as
  the starting values for those slots.  A numeric
  \code{start} argument of the appropriate length is used
  as the starting value of \code{theta}.}

  \item{optimizer}{which optimizer to use in the second
  phase of optimization.  The default is
  \code{\link{NelderMead}} and the alternative is
  \code{\link{bobyqa}}.  For difficult model fits we have
  found \code{\link{NelderMead}} to be more reliable but
  occasionally slower than \code{\link{bobyqa}}.}

  \item{mustart}{optional starting values on the scale of
  the conditional mean, as in \code{\link[stats]{glm}}; see
  there for details.}

  \item{etastart}{optional starting values on the scale of
  the unbounded predictor as in \code{\link[stats]{glm}};
  see there for details.}

  \item{tolPwrss}{numeric scalar - the tolerance for
  declaring convergence in the penalized iteratively
  weighted residual sum-of-squares step.  Defaults to
  1e-10.}

  \item{\dots}{other potential arguments.  A \code{method}
  argument was used in earlier versions of the package. Its
  functionality has been replaced by the \code{nAGQ}
  argument.}

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
  allowed during the optimization. Default is \code{900}}.
  \item{\code{tol}:}{a positive number specifying the
  convergence tolerance, currently only for the PWRSS
  iterations in \code{glmer()}.  Default is
  \code{0.000001}}.  }}

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
}
\value{
  An object of class \code{"\linkS4class{merMod}"}, for
  which many methods are available.  See there for details.
}
\description{
  Fit a generalized linear mixed model (GLMM)
}
\details{
  Fit a generalized linear mixed model, which incorporates
  both fixed-effects parameters and random effects in a
  linear predictor, via maximum likelihood. The linear
  predictor is related to the conditional mean of the
  response through the inverse link function defined in the
  GLM \code{family}.

  The expression for the likelihood of a mixed-effects
  model is an integral over the random effects space. For a
  linear mixed-effects model (LMM), as fit by
  \code{\link{lmer}}, this integral can be evaluated
  exactly.  For a GLMM the integral must be approximated.
  The most reliable approximation for GLMMs with a single
  grouping factor for the random effects is adaptive
  Gauss-Hermite quadrature.  The \code{nAGQ} argument
  controls the number of nodes in the quadrature formula.
  A model with a single, scalar random-effects term could
  reasonably use up to 25 quadrature points per scalar
  integral.

  With vector-valued random effects the complexity of the
  Gauss-Hermite quadrature formulas increases dramatically
  with the dimension.  For a 3-dimensional vector-valued
  random effect \code{nAGQ=5} requires 93 evaluations of
  the GLM deviance per evaluation of the approximate GLMM
  deviance.  For 20-dimensional evaluations of the GLM
  deviance per evaluation of the approximate GLMM deviance.

  The default approximation is the Laplace approximation,
  corresponding to \code{nAGQ=1}.
}
\examples{
## generalized linear mixed model
(gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
              data = cbpp, family = binomial))
## using nAGQ=0L only gets close to the optimum
(gm1a <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
               cbpp, binomial, nAGQ = 0L))
## using nAGQ=9L provides a better evaluation of the deviance
(gm1b <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
               cbpp, binomial, nAGQ = 9L))
## check with nAGQ=25L
(gm1c <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
               cbpp, binomial, nAGQ = 25L))

## GLMM with individual-level variability (accounting for overdispersion)
cbpp$obs <- 1:nrow(cbpp)
(gm2 <- glmer(cbind(incidence, size - incidence) ~ period +
    (1 | herd) +  (1|obs),
              family = binomial, data = cbpp))
anova(gm1,gm2)
}
\seealso{
  The \code{\linkS4class{merMod}} class,
  \code{\link[stats]{glm}}
}
\concept{
  GLMM
}
\keyword{models}

