\name{lmer}
\title{Fit Linear Mixed-Effects Models}
\alias{lmer}
\concept{ LMM }
\description{
  Fit a linear mixed-effects model (LMM) to data, via restricted maximum
  likelihood (REML) or maximum
  likelihood.
}
\usage{
lmer(formula, data = NULL, REML = TRUE, control = lmerControl(),
     start = NULL, verbose = 0L, subset, weights, na.action,
     offset, contrasts = NULL, devFunOnly = FALSE)
}
\arguments{
  \item{formula}{a two-sided linear formula object describing both the
    fixed-effects and random-effects part of the model, with the
    response on the left of a \code{~} operator and the terms, separated
    by \code{+} operators, on the right.  Random-effects terms are
    distinguished by vertical bars (\code{|}) separating expressions
    for design matrices from grouping factors. By default, non-scalar
    random effects (where the design matrix has more than one column,
    e.g. \code{(1+x|f)}) are fitted with \emph{unstructured} (general
    positive semidefinite) covariance matrices.
    \itemize{
      \item Two vertical bars
      (\code{||}) can be used to specify multiple uncorrelated random
      effects for the same grouping variable. %----------
      With default settings, the \code{||}-syntax \emph{works
	only for design matrices containing numeric (continuous) predictors};
      to fit models with independent categorical effects, use
      \code{diag(f|g)} or set \code{options(lme4.doublevert.default = "diag_special")} (see
      \code{\link{getDoublevertDefault}}).
      \item Tags preceding a random effect term specify covariance
      structure:
      \itemize{
	\item \code{us} (default: `us(f|g)` is equivalent to `(f|g)`):
	unstructured, positive semi-definite
	\item \code{diag}: diagonal (all correlations set to zero).
	Specify \code{diag(f|g, hom = TRUE)}
	to fit a homogeneous diagonal covariance matrix
	\item \code{cs}: compound symmetric (all pairwise correlations
	set identical). Specify \code{cs(f|g, hom = TRUE)} for homogeneous
	variances.
	\item \code{ar1}: autoregressive order 1. Note that AR1 models
	are homogeneous by default; specify \code{ar1(f|g, hom = FALSE)}
	for heterogeneous variances.
      } % itemize (cov struct types)
    } % itemize (cov struct info)
  } % formula

  \item{data}{an optional data frame containing the variables named in
    \code{formula}.  By default the variables are taken from the
    environment from which \code{lmer} is called. While \code{data} is
    optional, the package authors \emph{strongly} recommend its use,
    especially when later applying methods such as \code{update} and
    \code{drop1} to the fitted model (\emph{such methods are not
    guaranteed to work properly if \code{data} is omitted}). If
    \code{data} is omitted, variables will be taken from the environment
    of \code{formula} (if specified as a formula) or from the parent
    frame (if specified as a character vector).}

  \item{REML}{logical scalar - Should the estimates be chosen to
    optimize the REML criterion (as opposed to the log-likelihood)?}

  \item{control}{a list (of correct class, resulting from
    \code{\link{lmerControl}()} or \code{\link{glmerControl}()}
    respectively) containing control parameters, including the nonlinear
    optimizer to be used and parameters to be passed through to the
    nonlinear optimizer, see the \code{*lmerControl} documentation for
    details.}

  \item{start}{a numeric vector or a named \link{list} with one optional
    component named \code{par} or \code{theta}, giving starting values
    for covariance parameters.  Numeric \code{start} is equivalent to
    \code{list(par = start)}.  Parameters corresponding to unstructured
    covariance matrices are on the scale of the Cholesky factor of the
    relative covariance matrix.  By default, all relative covariance
    matrices are identity matrices.}

  \item{verbose}{integer scalar.  If \code{> 0} verbose output is
    generated during the optimization of the parameter estimates.  If
    \code{> 1} verbose output is generated during the individual
    penalized iteratively reweighted least squares (PIRLS) steps.}

  \item{subset}{an optional expression indicating the subset of the rows
    of \code{data} that should be used in the fit. This can be a logical
    vector, or a numeric vector indicating which observation numbers are
    to be included, or a character vector of the row names to be
    included.  All observations are included by default.}

  \item{weights}{an optional vector of \sQuote{prior weights} to be used
    in the fitting process.  Should be \code{NULL} or a numeric vector.
    Prior \code{weights} are \emph{not} normalized or standardized in
    any way.  In particular, the diagonal of the residual covariance
    matrix is the squared residual standard deviation parameter
    \code{\link{sigma}} times the vector of inverse \code{weights}.
    Therefore, if the \code{weights} have relatively large magnitudes,
    then in order to compensate, the \code{\link{sigma}} parameter will
    also need to have a relatively large magnitude.}

  \item{na.action}{a function that indicates what should happen when the
    data contain \code{NA}s.  The default action (\code{na.omit},
    inherited from the 'factory fresh' value of
    \code{getOption("na.action")}) strips any observations with any
    missing values in any variables.}

  \item{offset}{this can be used to specify an \emph{a priori} known
    component to be included in the linear predictor during
    fitting. This should be \code{NULL} or a numeric vector of length
    equal to the number of cases.  One or more \code{\link{offset}}
    terms can be included in the formula instead or as well, and if more
    than one is specified their sum is used.  See
    \code{\link{model.offset}}.}

  \item{contrasts}{an optional list. See the \code{contrasts.arg} of
    \code{model.matrix.default}.}

  \item{devFunOnly}{logical - return only the deviance evaluation
    function. Note that because the deviance function operates on
    variables stored in its environment, it may not return
    \emph{exactly} the same values on subsequent calls (but the results
    should always be within machine tolerance).}
}
\value{
  An object of class \code{\link[=merMod-class]{merMod}} (more specifically,
  an object of \emph{subclass} \code{lmerMod}), for which many methods
  are available (e.g. \code{methods(class="merMod")})
}
\note{
  In earlier version of the \pkg{lme4} package, a \code{method} argument was
  used.  Its functionality has been replaced by the \code{REML} argument.

  Also, \code{lmer(.)} allowed a \code{family} argument (to effectively
  switch to \code{glmer(.)}).  This has been deprecated in summer 2013,
  and been disabled in spring 2019.
}
\details{
  \itemize{
    \item If the \code{formula} argument is specified as a character
      vector, the function will attempt to coerce it to a formula.
      However, this is not recommended (users who want to construct
      formulas by pasting together components are advised to use
      \code{\link{as.formula}} or \code{\link{reformulate}}); model fits
      will work but subsequent methods such as \code{\link{drop1}},
      \code{\link{update}}, etc. may fail.

    \item When handling perfectly collinear predictor variables
      (i.e. fixed-effect design matrices of less than full rank),
      \code{[gn]lmer} is not as sophisticated
      as modeling frameworks such as
      \code{\link{lm}} and \code{\link{glm}}. While it does
      automatically drop collinear variables (with a message
      rather than a warning), it does not automatically fill
      in \code{NA} values for the dropped coefficients;
      these can be added via \code{fixef(fitted.model, add.dropped=TRUE)}.
      This information can also be retrieved via
      \code{attr(getME(fitted.model, "X"), "col.dropped")}.
    
    \item the deviance function returned when \code{devFunOnly} is
      \code{TRUE} takes a single numeric vector argument which defines
      the scaled variance-covariance matrices of the random effects.
      \itemize{
	\item In the case of unstructured covariances, this vector
	is directly mapped to the \code{theta} vector, which represents
	the unique non-zero values in the Cholesky factor of the
	(scaled) covariance matrix. For models with only simple
	(intercept-only) random effects, \code{par} (and thus
	\code{theta}) is a vector of the
	standard deviations of the random effects.  For more complex or
	multiple random effects, running \code{getME(.,"par")} or
	(equivalently) \code{getME(., "theta")} to retrieve the
	\code{theta} 
	vector for a fitted model and examining
	the names of the vector is probably the easiest way to determine
	the correspondence between the elements of the \code{theta} vector
	and elements of the lower triangles of the Cholesky factors of the
	random effects.
	\item For structured covariances, the \code{getTheta} method
	translates the parameter vector to the \code{theta}
	(Cholesky-factor element) scale for internal use. The parameter
	vector is usually composed of a set of standard-deviation values
	(one if \code{hom = TRUE} or many if \code{hom = FALSE}),
	followed by one or more parameters that determine the correlation matrix.
      } % itemize (covariance parameterizations)	
    } % itemize (details)
} % details
\seealso{
  \code{\link[stats]{lm}} for linear models;
  \code{\link{glmer}} for generalized linear; and
  \code{\link{nlmer}} for nonlinear mixed models.
  
  \code{\link{plot.merMod}} for plot diagnostics.
}
\examples{
## linear mixed models - reference values from older code
(fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy))
summary(fm1) # (with its own print method; see class?merMod % ./merMod-class.Rd
plot(fm1) # plotting the model diagnostics; see ?plot.merMod

str(terms(fm1))
stopifnot(identical(terms(fm1, fixed.only=FALSE),
                    terms(model.frame(fm1))))
attr(terms(fm1, FALSE), "dataClasses") # fixed.only=FALSE needed for dataCl.

## Maximum Likelihood (ML), and "monitor" iterations via 'verbose':
fm1_ML <- update(fm1, REML=FALSE, verbose = 1)
(fm2 <- lmer(Reaction ~ Days + (Days || Subject), sleepstudy))
anova(fm1, fm2)
sm2 <- summary(fm2)
print(fm2, digits=7, ranef.comp="Var") # the print.merMod()         method
print(sm2, digits=3, corr=FALSE)       # the print.summary.merMod() method

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
