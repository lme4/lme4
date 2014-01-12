\docType{class}
\name{merMod-class}
\title{Class "merMod" of Fitted Mixed-Effect Models}
\alias{anova.merMod}
\alias{coef.merMod}
\alias{deviance.merMod}
\alias{fitted.merMod}
\alias{formula.merMod}
\alias{glmerMod-class}
\alias{lmerMod-class}
\alias{logLik.merMod}
\alias{merMod}
\alias{merMod-class}
\alias{model.frame.merMod}
\alias{model.matrix.merMod}
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
\S3method{terms}{merMod}(x, fixed.only = TRUE, \dots)
\S3method{vcov}{merMod}(object, correlation = TRUE, sigm = sigma(object),
use.hessian = NULL, \dots)

\S3method{print}{merMod}(x, digits = max(3, getOption("digits") - 3),
      correlation = NULL, symbolic.cor = FALSE,
      signif.stars = getOption("show.signif.stars"), ranef.comp = "Std.Dev.", ...)

\S3method{summary}{merMod}(object, correlation = , use.hessian = NULL, \dots)
\S3method{print}{summary.merMod}(x, digits = max(3, getOption("digits") - 3),
      correlation = NULL, symbolic.cor = FALSE,
      signif.stars = getOption("show.signif.stars"),
      ranef.comp = c("Variance", "Std.Dev."), show.resids = TRUE, ...)
}
\arguments{
  \item{object}{an \R object of class \code{\linkS4class{merMod}}, i.e.,
    as resulting from \code{\link{lmer}()}, or \code{\link{glmer}()},
    etc.}
  \item{x}{an \R object of class \code{merMod} or \code{summary.merMod},
    respectively, the latter resulting from \code{summary(<merMod>)}.}
  \item{refit}{logical indicating if objects of class \code{lmerMod} should be
    refitted with ML before producing the anova table, default is \code{TRUE}.}
  \item{model.names}{character vectors of model names to be used in the anova table.}
  \item{fixed.only}{logical indicating if only the fixed effects
    \code{\link{terms}} are sought, defaults to true.  If false, all
    terms, including random ones are returned.}
  \item{correlation}{(logical) for \code{vcov}, indicates whether the
    correlation matrix as well as the variance-covariance matrix is
    desired; for \code{print.summary.merMod}, indicates whether the
    correlation matrix of the fixed-effects parameters should be
    printed.}
  \item{use.hessian}{(logical) indicates whether to use the
    finite-difference Hessian of the deviance function to compute
    standard errors of the fixed effects, rather estimating
    based on internal information about the inverse of
    the model matrix (see \code{\link{getME}(.,"RX")}).
    The default is to to use the Hessian whenever the
    fixed effect parameters are arguments to the deviance
    function, rather than being profiled out (i.e. for
    GLMMs with \code{nAGQ>0}).

    \code{use.hessian=FALSE} is backward-compatible with older versions
    of \code{lme4}, but may give less accurate SE estimates when the
    estimates of the fixed-effect and random-effect parameters
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
  \item{show.resids}{should the quantiles of the scaled residuals be printed?}
  \item{\dots}{potentially further arguments passed from other methods.}
}
\section{Objects from the Class}{
  Objects of class \code{merMod} are created by calls to
  \code{\link{lmer}}, \code{\link{glmer}} or \code{\link{nlmer}}.
}
\section{S3 methods}{
The following S3 methods with arguments given above exist (this list is currently not complete):
  \describe{
  \item{\code{summary}:}{Computes and returns a list of summary statistics of the
      fitted model, the amount of output can be controlled via the \code{print} method,
      see also \code{\link{summary}}.}
   \item{\code{print.summary}:}{Controls the output for the summary method.}
  \item{\code{vcov}:}{Calculate variance-covariance matrix of the \emph{fixed}
     effect terms, see also \code{\link{vcov}}.}
  \item{\code{anova}:}{returns the sequential decomposition of the contributions of
      fixed-effects terms or, for multiple arguments, model comparison statistics.
      For objects of class \code{lmerMod} the default behavior is to refit the models
      with LM if fitted with \code{REML = TRUE}, this can be controlled via the
      \code{refit} argument. See also \code{\link{anova}}.}
      }
}
\seealso{
  \code{\link{lmer}}, \code{\link{glmer}},
  \code{\link{nlmer}}, \code{\linkS4class{merPredD}},
  \code{\linkS4class{lmerResp}},
  \code{\linkS4class{glmResp}},
  \code{\linkS4class{nlsResp}}
}
\examples{
showClass("merMod")
methods(class="merMod")## over 30  (S3) methods available

## -> example(lmer)  for an example of vcov.merMod()
}
\keyword{classes}
