\docType{class}
\name{merMod-class}
\title{Class "merMod" of Fitted Mixed-Effect Models}
\alias{anova.merMod}
\alias{as.function.merMod}
\alias{coef.merMod}
\alias{deviance.merMod}
\alias{df.residual.merMod}
\alias{extractAIC.merMod}
\alias{family.merMod}
\alias{fitted.merMod}
\alias{formula.merMod}
\alias{glmerMod-class}
\alias{lmerMod-class}
\alias{logLik.merMod}
\alias{merMod}
\alias{merMod-class}
\alias{model.frame.merMod}
\alias{model.matrix.merMod}
\alias{ngrps.merMod}
\alias{nobs.merMod}
\alias{nobs}
\alias{nlmerMod-class}
\alias{print.merMod}
\alias{print.summary.merMod}
\alias{show,merMod-method}
\alias{show.merMod}
\alias{show.summary.merMod}
\alias{summary.merMod}
\alias{summary.summary.merMod}
\alias{terms.merMod}
\alias{update.merMod}
\alias{vcov.merMod}
\alias{vcov.summary.merMod}
\alias{weights.merMod}
\alias{REMLcrit}
\description{
  A mixed-effects model is represented as a
  \code{\linkS4class{merPredD}} object and a response
  module of a class that inherits from class
  \code{\linkS4class{lmResp}}.  A model with a
  \code{\linkS4class{lmerResp}} response has class \code{lmerMod}; a
  \code{\linkS4class{glmResp}} response has class \code{glmerMod}; and a
  \code{\linkS4class{nlsResp}} response has class \code{nlmerMod}.
}
\usage{
\S3method{anova}{merMod}(object, ..., refit = TRUE, model.names=NULL)
\S3method{as.function}{merMod}(x, ...)
\S3method{coef}{merMod}(object, ...)
\S3method{deviance}{merMod}(object, REML = NULL, ...)
REMLcrit(object)
\S3method{extractAIC}{merMod}(fit, scale = 0, k = 2, ...)
\S3method{family}{merMod}(object, ...)
\S3method{formula}{merMod}(x, fixed.only = FALSE, random.only = FALSE, ...)
\S3method{fitted}{merMod}(object, ...)
\S3method{logLik}{merMod}(object, REML = NULL, ...)
\S3method{nobs}{merMod}(object, ...)
\S3method{ngrps}{merMod}(object, ...)
\S3method{terms}{merMod}(x, fixed.only = TRUE, random.only = FALSE, \dots)
\S3method{vcov}{merMod}(object, correlation = TRUE, sigm = sigma(object),
     use.hessian = NULL, \dots)
\S3method{model.frame}{merMod}(formula, fixed.only = FALSE, ...)
\S3method{model.matrix}{merMod}(object, type = c("fixed", "random", "randomListRaw"), ...)
\S3method{print}{merMod}(x, digits = max(3, getOption("digits") - 3),
      correlation = NULL, symbolic.cor = FALSE,
      signif.stars = getOption("show.signif.stars"), ranef.comp = "Std.Dev.", ...)

\S3method{summary}{merMod}(object, correlation = , use.hessian = NULL, \dots)
\S3method{print}{summary.merMod}(x, digits = max(3, getOption("digits") - 3),
      correlation = NULL, symbolic.cor = FALSE,
      signif.stars = getOption("show.signif.stars"),
      ranef.comp = c("Variance", "Std.Dev."), show.resids = TRUE, ...)
\S3method{update}{merMod}(object, formula., ..., evaluate = TRUE)
\S3method{weights}{merMod}(object, type = c("prior", "working"), ...)
}
\arguments{
  \item{object}{an \R object of class \code{\linkS4class{merMod}}, i.e.,
    as resulting from \code{\link{lmer}()}, or \code{\link{glmer}()},
    etc.}
  \item{x}{an \R object of class \code{merMod} or \code{summary.merMod},
    respectively, the latter resulting from \code{summary(<merMod>)}.}
  \item{fit}{an \R object of class \code{\linkS4class{merMod}}.}
  \item{formula}{in the case of \code{model.frame}, a
    \code{\linkS4class{merMod}} object.}
  \item{refit}{logical indicating if objects of class \code{lmerMod} should be
    refitted with ML before comparing models. The default is
    \code{TRUE} to prevent the common mistake of inappropriately
    comparing REML-fitted models with different fixed effects,
    whose likelihoods are not directly comparable.}
  \item{model.names}{character vectors of model names to be used in the
    anova table.}
  \item{scale}{Not currently used (see \code{\link{extractAIC}}).}
  \item{k}{see \code{\link{extractAIC}}.}
  \item{REML}{Logical. If \code{TRUE}, return the restricted log-likelihood
    rather than the log-likelihood.  If \code{NULL} (the default),
    set \code{REML} to \code{isREML(object)} (see \code{\link{isREML}}).}
  \item{fixed.only}{logical indicating if only the fixed effects
    components (terms or formula elements) are sought.  If false, all
    components, including random ones, are returned.}
  \item{random.only}{complement of \code{fixed.only}; indicates
    whether random components only are sought. (Trying to specify
    \code{fixed.only} and \code{random.only} at the same time
    will produce an error.)}
  \item{correlation}{(logical)
    for \code{vcov}, indicates whether the correlation matrix as well as
    the variance-covariance matrix is desired;
    for \code{summary.merMod}, indicates whether the correlation matrix
    should be computed and stored along with the covariance;
    for \code{print.summary.merMod}, indicates whether the correlation
    matrix of the fixed-effects parameters should be printed.  In the
    latter case, when \code{NULL} (the default), the correlation matrix
    is printed when it has been computed by \code{summary(.)}, and when
    \eqn{p <= 20}.}% and '20' can be changed by  options(lme4.summary.cor.max = <n>)
  \item{use.hessian}{(logical) indicates whether to use the
    finite-difference Hessian of the deviance function to compute
    standard errors of the fixed effects, rather estimating
    based on internal information about the inverse of
    the model matrix (see \code{\link{getME}(.,"RX")}).
    The default is to to use the Hessian whenever the
    fixed effect parameters are arguments to the deviance
    function (i.e. for GLMMs with \code{nAGQ>0}), and to use
    \code{\link{getME}(.,"RX")} whenever the fixed effect parameters are
    profiled out (i.e. for GLMMs with \code{nAGQ==0} or LMMs).

    \code{use.hessian=FALSE} is backward-compatible with older versions
    of \code{lme4}, but may give less accurate SE estimates when the
    estimates of the fixed-effect (see \code{\link{getME}(.,"beta")})
    and random-effect (see \code{\link{getME}(.,"theta")}) parameters
    are correlated.
    }
  \item{sigm}{the residual standard error; by default \code{\link{sigma}(object)}.}
  \item{digits}{number of significant digits for printing}
  \item{symbolic.cor}{should a symbolic encoding of the fixed-effects correlation
    matrix be printed? If so, the \code{\link{symnum}} function is used.}
  \item{signif.stars}{(logical) should significance stars be used?}
  \item{ranef.comp}{character vector of length one or two, indicating
    if random-effects parameters should be reported on the variance and/or
    standard deviation scale.}
  \item{show.resids}{should the quantiles of the scaled residuals be
    printed?}
  \item{formula.}{see \code{\link{update.formula}}.}
  \item{evaluate}{see \code{\link{update}}.}
  \item{type}{For \code{weights}, type of weights to be returned; either \code{"prior"} for
    the initially supplied weights or \code{"working"} for the weights
    at the final iteration of the penalized iteratively reweighted least
    squares algorithm. For \code{model.matrix}, type of model matrix to
    return (one of \code{fixed} giving the fixed effects model matrix,
    \code{random} giving the random effects model matrix, or
    \code{randomListRaw} giving a list of the raw random effects model
    matrices associated with each random effects term).}
  \item{\dots}{potentially further arguments passed from other methods.}
}
\section{Objects from the Class}{
  Objects of class \code{merMod} are created by calls to
  \code{\link{lmer}}, \code{\link{glmer}} or \code{\link{nlmer}}.
}
\section{S3 methods}{
The following S3 methods with arguments given above exist (this list is currently not complete):
\describe{
  %% TODO:  document differences between update and update.merMod
  \item{\code{anova}:}{returns the sequential decomposition of the contributions of
    fixed-effects terms or, for multiple arguments, model comparison statistics.
    For objects of class \code{lmerMod} the default behavior is to refit the models
    with ML if fitted with \code{REML = TRUE}, this can be controlled via the
    \code{refit} argument. See also \code{\link{anova}}.}
  \item{\code{as.function}:}{returns the deviance function, the same as
    \code{\link{lmer}(*, devFunOnly=TRUE)}, and \code{\link{mkLmerDevfun}()}
    or \code{\link{mkGlmerDevfun}()}, respectively.}
  \item{\code{coef}:}{Computes the sum of the random and fixed effects
    coefficients for each explanatory variable for each level of each
    grouping factor.}
  \item{\code{extractAIC}:}{Computes the (generalized) Akaike An
    Information Criterion. If \code{isREML(fit)}, then \code{fit} is
    refitted using maximum likelihood.}
  \item{\code{family}:}{\code{\link{family}} of fitted
    GLMM. (\emph{Warning:} this accessor may not work properly with
    customized families/link functions.)}
  \item{\code{fitted}:}{Fitted values, given the conditional modes of
    the random effects.  For more flexible access to fitted values, use
    \code{\link{predict.merMod}}.}
  \item{\code{logLik}:}{Log-likelihood at the fitted value of the
    parameters.  Note that for GLMMs, the returned value is only
    proportional to the log probability density (or distribution) of the
    response variable.  See \code{\link{logLik}}. }
  \item{\code{model.frame}:}{returns the \code{frame} slot of \code{\linkS4class{merMod}}.}
  \item{\code{model.matrix}:}{returns the fixed effects model matrix.}
  \item{\code{nobs}, \code{ngrps}:}{Number of observations and vector of
    the numbers of levels in each grouping factor.  See \code{\link{ngrps}}.}
  \item{\code{summary}:}{Computes and returns a list of summary statistics of the
      fitted model, the amount of output can be controlled via the \code{print} method,
      see also \code{\link{summary}}.}
   \item{\code{print.summary}:}{Controls the output for the summary
     method.}
  \item{\code{vcov}:}{Calculate variance-covariance matrix of the \emph{fixed}
     effect terms, see also \code{\link{vcov}}.}
  \item{\code{update}:}{See \code{\link{update}}.}
      }
}
\section{Deviance and log-likelihood of GLMMs}{

  One must be careful when defining the deviance of a GLM.  For example,
  should the deviance be defined as minus twice the log-likelihood or
  does it involve subtracting the deviance for a saturated model?  To
  distinguish these two possibilities we refer to absolute deviance
  (minus twice the log-likelihood) and relative deviance (relative to a
  saturated model, e.g. Section 2.3.1 in McCullagh and Nelder 1989).
  With GLMMs however, there is an additional complication involving the
  distinction between the likelihood and the conditional likelihood.
  The latter is the likelihood obtained by conditioning on the estimates
  of the conditional modes of the spherical random effects coefficients,
  whereas the likelihood itself (i.e. the unconditional likelihood)
  involves integrating out these coefficients. The following table
  summarizes how to extract the various types of deviance for a
  \code{glmerMod} object:

  \tabular{lrr}{
             \tab conditional              \tab unconditional \cr
    relative \tab \code{deviance(object)}  \tab NA in \code{lme4} \cr
    absolute \tab \code{object@resp$aic()} \tab \code{-2*logLik(object)}
  }

  This table requires two caveats:

  \itemize{

    \item If the link function involves a scale parameter
    (e.g. \code{Gamma}) then \code{object@resp$aic() - 2 * getME(object,
    "devcomp")$dims["useSc"]} is required for the absolute-conditional
    case.
    \item If adaptive Gauss-Hermite quadrature is used, then
    \code{logLik(object)} is currently only proportional to the
    absolute-unconditional log-likelihood.

  }

  For more information about this topic see the \code{misc/logLikGLMM}
  directory in the package source.

  %% We define the conditional likelihood as the density of the response
  %% variable, \eqn{y}, conditional on the spherical random effects,
  %% \eqn{u}:

  %% \deqn{f_{\theta,\beta}(y|u)}{f(y|u)}

  %% which depends on the fixed effect and covariance parameters,
  %% \eqn{\beta} and \eqn{\theta}.  We define the likelihood itself as:

  %% \deqn{f_{\theta,\beta}(y) = \int f_{\theta,\beta}(y|u) f(u)}{f(y) = integral(f(y|u)f(u))}

  %% where \eqn{f(u)} is an independent multivariate normal.  It is not
  %% typically possible to evaluate this integral for GLMMs.  In
  %% \code{lme4} we approximate it using the Laplace approximation and
  %% adaptive Gauss-Hermite quadrature.

  %% For canonical links, the Laplace approximation to the log-likelihood
  %% is:

  %% \deqn{log f_{\theta,\beta}(y | u)-0.5|u|^2 - log|L_{\theta}|}
  %% {log f(y|u) - 0.5|u|^2 - log|L|}

  %% where \eqn{|u|^2} is the squared length of the vector of spherical
  %% random effects, \eqn{u}, and \eqn{|L_{\theta}|}{|L|} is the
  %% determinant of the Cholesky factor for the cross-product matrix for
  %% the random effects (\code{\link{getME}(., "L")}).  Note that in the
  %% non-canonical link case, the above equation is technically an
  %% approximation to the Laplace approximation.

  %% For models with a single scalar random effect one may also obtain an
  %% approximate log-likelihood using adaptive Gauss-Hermite quadrature as:

  %% \deqn{-\log|\bm L| - \frac{q}{2}\log(2\pi) + \sum_j^q \log \left[\sum_k^{N_{AGQ}} w_k f(y, u_{jk}) e^{-\frac{1}{2}u_{jk}^2 - l_k}\right]}

  %% where \eqn{N_{AGQ}} is the number of quadrature nodes, \eqn{q} is the
  %% number of levels of the random effect, the \eqn{w_k} and \eqn{l_k} are
  %% weights and log of the normal density associated with the \eqn{k}th
  %% quadrature node (\code{\link{GHrule}}), and \eqn{u_{jk}} is the value
  %% of the spherical random effect coefficient associated with the
  %% \eqn{j}th level of the grouping factor and \eqn{k}th quadrature node.

}
\section{Slots}{
  \describe{
    \item{\code{resp}:}{A reference class object for an \pkg{lme4}
      response module (\code{\link{lmResp-class}}).}
    \item{\code{Gp}:}{See \code{\link{getME}}.}
    \item{\code{call}:}{The matched call.}
    \item{\code{frame}:}{The model frame containing all of the variables
      required to parse the model formula.}
    \item{\code{flist}:}{See \code{\link{getME}}.}
    \item{\code{cnms}:}{See \code{\link{getME}}.}
    \item{\code{lower}:}{See \code{\link{getME}}.}
    \item{\code{theta}:}{Covariance parameter vector.}
    \item{\code{beta}:}{Fixed effects coefficients.}
    \item{\code{u}:}{Conditional model of spherical random effects
      coefficients.}
    \item{\code{devcomp}:}{See \code{\link{getME}}.}
    \item{\code{pp}:}{A reference class object for an \pkg{lme4}
      predictor module (\code{\link{merPredD-class}}).}
    \item{\code{optinfo}:}{List containing information about the
      nonlinear optimization.}
  }
}
\seealso{
  \code{\link{lmer}}, \code{\link{glmer}},
  \code{\link{nlmer}}, \code{\linkS4class{merPredD}},
  \code{\linkS4class{lmerResp}},
  \code{\linkS4class{glmResp}},
  \code{\linkS4class{nlsResp}}

  Other methods for \code{merMod} objects documented elsewhere include:
  \code{\link{fortify.merMod}}, \code{\link{drop1.merMod}},
  \code{\link{isLMM.merMod}}, \code{\link{isGLMM.merMod}},
  \code{\link{isNLMM.merMod}}, \code{\link{isREML.merMod}},
  \code{\link{plot.merMod}}, \code{\link{predict.merMod}},
  \code{\link{profile.merMod}}, \code{\link{ranef.merMod}},
  \code{\link{refit.merMod}}, \code{\link{refitML.merMod}},
  \code{\link{residuals.merMod}}, \code{\link{sigma.merMod}},
  \code{\link{simulate.merMod}}, \code{\link{summary.merMod}}.
}
\examples{
showClass("merMod")
methods(class="merMod")## over 30  (S3) methods available

## -> example(lmer)  for an example of vcov.merMod()
}
\keyword{classes}
