\name{glmer}
\alias{glmer}
\title{Fit Generalized Linear Mixed-Effects Models}
\usage{
  glmer(formula, data = NULL, family = gaussian,
    control = glmerControl(), start = NULL, verbose = 0L,
    nAGQ = 1L, subset, weights, na.action, offset,
    contrasts = NULL, mustart, etastart,
    devFunOnly = FALSE, ...)
}
\arguments{
  \item{family}{a GLM family, see \code{\link[stats]{glm}}
  and \code{\link[stats]{family}}.}

  \item{nAGQ}{integer scalar - the number of points per
  axis for evaluating the adaptive Gauss-Hermite
  approximation to the log-likelihood.  Defaults to 1,
  corresponding to the Laplace approximation.  Values
  greater than 1 produce greater accuracy in the evaluation
  of the log-likelihood at the expense of speed.  A value
  of zero uses a faster but less exact form of parameter
  estimation for GLMMs by optimizing the random effects and
  the fixed-effects coefficients in the penalized
  iteratively reweighted least squares step.}

  \item{start}{a named list of starting values for the
  parameters in the model, or a numeric vector. A numeric
  \code{start} argument will be used as the starting value
  of \code{theta}.  If \code{start} is a list, the
  \code{theta} element (a numeric vector) is used as the
  starting value for the first optimization step (default=1
  for diagonal elements and 0 for off-diagonal elements of
  the lower Cholesky factor); the fitted value of
  \code{theta} from the first step, plus
  \code{start[["fixef"]]}, are used as starting values for
  the second optimization step.  If \code{start} has both
  \code{fixef} and \code{theta} elements, the first
  optimization step is skipped. For more details or finer
  control of optimization, see \code{\link{modular}}.}

  \item{mustart}{optional starting values on the scale of
  the conditional mean, as in \code{\link[stats]{glm}}; see
  there for details.}

  \item{etastart}{optional starting values on the scale of
  the unbounded predictor as in \code{\link[stats]{glm}};
  see there for details.}

  \item{\dots}{other potential arguments.  A \code{method}
  argument was used in earlier versions of the package. Its
  functionality has been replaced by the \code{nAGQ}
  argument.}

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

  \item{control}{a list (of correct class, resulting from
  \code{\link{lmerControl}()} or
  \code{\link{glmerControl}()} respectively) containing
  control parameters, including the nonlinear optimizer to
  be used and parameters to be passed through to the
  nonlinear optimizer, see the \code{*lmerControl}
  documentation for details.}

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
}
\value{
  An object of class \code{glmerMod}, for which many
  methods are available (e.g.
  \code{methods(class="glmerMod")})
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
library(lattice)
xyplot(incidence/size ~ period|herd, cbpp, type=c('g','p','l'),
       layout=c(3,5), index.cond = function(x,y)max(y))
(gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
              data = cbpp, family = binomial))
## using nAGQ=0 only gets close to the optimum
(gm1a <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
               cbpp, binomial, nAGQ = 0))
## using  nAGQ = 9  provides a better evaluation of the deviance
## Currently the internal calculations use the sum of deviance residuals,
## which is not directly comparable with the nAGQ=0 or nAGQ=1 result.
(gm1a <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
               cbpp, binomial, nAGQ = 9))

## GLMM with individual-level variability (accounting for overdispersion)
## For this data set the model is the same as one allowing for a period:herd
## interaction, which the plot indicates could be needed.
cbpp$obs <- 1:nrow(cbpp)
(gm2 <- glmer(cbind(incidence, size - incidence) ~ period +
    (1 | herd) +  (1|obs),
              family = binomial, data = cbpp))
anova(gm1,gm2)

## glmer and glm log-likelihoods are consistent
gm1Devfun <- update(gm1,devFunOnly=TRUE)
gm0 <- glm(cbind(incidence, size - incidence) ~ period,
           family = binomial, data = cbpp)
## evaluate GLMM deviance at RE variance=theta=0, beta=(GLM coeffs)
gm1Dev0 <- gm1Devfun(c(0,coef(gm0)))
## compare
stopifnot(all.equal(gm1Dev0,c(-2*logLik(gm0))))
}
\seealso{
  \code{\link{lmer}} (for details on formulas and
  parameterization); \code{\link[stats]{glm}}
}
\concept{
  GLMM
}
\keyword{models}

