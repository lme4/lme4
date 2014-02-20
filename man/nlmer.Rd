\name{nlmer}
\alias{nlmer}
\title{Fit Nonlinear Mixed-Effects Models}
\usage{
  nlmer(formula, data = NULL, control = nlmerControl(),
    start = NULL, verbose = 0L, nAGQ = 1L, subset, weights,
    na.action, offset, contrasts = NULL,
    devFunOnly = FALSE, ...)
}
\arguments{
  \item{formula}{a three-part \dQuote{nonlinear mixed model} formula,
    of the form \code{resp ~ Nonlin(...) ~ fixed + random},
    where the third part is similar to the RHS formula of, e.g.,
    \code{\link{lmer}}. %% (see detailed documentation)  FIXME!  where would that be?
    
    Currently, the \code{Nonlin(..)} formula part must not only return a
    numeric vector, but also must have a \code{"gradient"} attribute, a
    \code{\link{matrix}}.  The \code{\link{SSbiexp}},
    \code{\link{SSlogis}}, etc, see \code{\link{selfStart}} functions
    provide this (and more).  Alternatively, you can use
    \code{\link{deriv}()} to automatically produce such functions or
    expressions.}

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
\description{
  Fit a nonlinear mixed-effects model
}
\details{
  Fit nonlinear mixed-effects models, such as those used in
  population pharmacokinetics.
}
\note{
  Adaptive Gauss-Hermite quadrature (\code{nAGQ}>1) is not
  currently implemented for \code{nlmer}. Several other
  methods, such as simulation or prediction with new data,
  are unimplemented or very lightly tested.
}
\examples{
## nonlinear mixed models --- 3-part formulas ---
## 1. basic nonlinear fit. Use stats::SSlogis for its
## implementation of the 3-parameter logistic curve.
## "SS" stands for "self-starting logistic", but the
## "self-starting" part is not currently used by nlmer ... 'start' is
## necessary
startvec <- c(Asym = 200, xmid = 725, scal = 350)
(nm1 <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree,
             Orange, start = startvec))
## 2. re-run with "quick and dirty" PIRLS step
(nm1a <- update(nm1, nAGQ = 0L))
## 3. Fit the same model with a user-built function:
## a. Define formula
nform <- ~Asym/(1+exp((xmid-input)/scal))
## b. Use deriv() to construct function:
nfun <- deriv(nform,namevec=c("Asym","xmid","scal"),
              function.arg=c("input","Asym","xmid","scal"))
nm1b <- update(nm1,circumference ~ nfun(age, Asym, xmid, scal)  ~ Asym | Tree)
## 4. User-built function without using derivs():
##    derivatives could be computed more efficiently
##    by pre-computing components, but these are essentially
##    the gradients as one would derive them by hand
nfun2 <- function(input, Asym, xmid, scal) {
    value <- Asym/(1+exp((xmid-input)/scal))
    grad <- cbind(Asym=1/(1+exp((xmid-input)/scal)),
              xmid=-Asym/(1+exp((xmid-input)/scal))^2*1/scal*
                    exp((xmid-input)/scal),
              scal=-Asym/(1+exp((xmid-input)/scal))^2*
                     -(xmid-input)/scal^2*exp((xmid-input)/scal))
    attr(value,"gradient") <- grad
    value
}
stopifnot(all.equal(attr(nfun(2,1,3,4),"gradient"),
                    attr(nfun(2,1,3,4),"gradient")))
nm1c <- update(nm1,circumference ~ nfun2(age, Asym, xmid, scal)  ~ Asym | Tree)
}
\keyword{models}

