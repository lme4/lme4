\name{merMod-class}
\title{Class "merMod" of Fitted Mixed-Effect Models}
\Rdversion{1.1}
\docType{class}
\alias{merMod-class}
\alias{anova,merMod-method}
\alias{coef,merMod-method}
%\alias{devcomp,merMod-method}
\alias{deviance,merMod-method}
\alias{fitted,merMod-method}
%\alias{fixef,merMod-method}
\alias{formula,merMod-method}
%\alias{getL,merMod-method}
\alias{logLik,merMod-method}
\alias{model.frame,merMod-method}
\alias{model.matrix,merMod-method}
\alias{print,merMod-method}
%\alias{profile,merMod-method}
%\alias{ranef,merMod-method}
\alias{show,merMod-method}
%\alias{sigma,merMod-method}
\alias{simulate,merMod-method}
\alias{summary,merMod-method}
\alias{terms,merMod-method}
\alias{update,merMod-method}
\alias{VarCorr,merMod-method}
\alias{vcov,merMod-method}
% FIXME ? : "summary.mer" is not documented
\alias{print,summary.mer-method}
\alias{show,summary.mer-method}
\alias{summary,summary.mer-method}
\alias{vcov,summary.mer-method}
%
\description{
  A mixed-effects model represented as a \code{\linkS4class{reModule}}
  for the random effects, an \code{\linkS4class{feModule}} for the
  fixed-effects parameters and an \code{\linkS4class{respModule}} module
  for the response.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("merMod", ...)}.
  More commonly the objects are created by calls to \code{\link{lmer}},
  \code{\link{glmer}} or \code{\link{nlmer}}.
}
\section{Slots}{
  \describe{
    \item{\code{call}:}{The matched call to the function that generated
      the object.}
    \item{\code{devcomp}:}{A list with two components:
      \dQuote{\code{cmp}}, a numeric vector of components used to
      evaluate the deviance, and \dQuote{\code{dims}}, an integer vector
      of dimensions and flags.}
    \item{\code{frame}:}{The model frame used to evaluate the model matrices.}
    \item{\code{re}:}{An \code{\linkS4class{reModule}} object
      representing the random effects values and various matrices
      associated with them.}
    \item{\code{fe}:}{An \code{\linkS4class{feModule}} object
      representing the fixed-effects parameters and various matrices
      associated with them.}
    \item{\code{resp}:}{An \code{\linkS4class{respModule}} object containing
	information on the response, the fitted values, the residuals
	and various forms of weights.  According to the class of this
	object the model is a linear mixed model (\code{\link{lmer}}), a
	generalized linear mixed model (\code{\link{glmer}}) or a
	nonlinear mixed model (\code{\link{nlmer}}).}
  }
}
\section{Methods}{
  \describe{
    \item{anova}{\code{signature(object = "merMod")}: ... }
    \item{coef}{\code{signature(object = "merMod")}: ... }
    \item{devcomp}{\code{signature(x = "merMod")}: ... }
    \item{deviance}{\code{signature(object = "merMod")}: ... }
    \item{fixef}{\code{signature(object = "merMod")}: ... }
    \item{formula}{\code{signature(x = "merMod")}: ... }
    \item{getL}{\code{signature(x = "merMod")}: ... }
    \item{logLik}{\code{signature(object = "merMod")}: ... }
    \item{model.frame}{\code{signature(formula = "merMod")}: ... }
    \item{model.matrix}{\code{signature(object = "merMod")}: ... }
    \item{print}{\code{signature(x = "merMod")}: ... }
    \item{profile}{\code{signature(fitted = "merMod")}: ... }
    \item{ranef}{\code{signature(object = "merMod")}: ... }
    \item{show}{\code{signature(object = "merMod")}: ... }
    \item{sigma}{\code{signature(object = "merMod")}: ... }
    \item{simulate}{\code{signature(object = "merMod")}: ... }
    \item{summary}{\code{signature(object = "merMod")}: ... }
    \item{terms}{\code{signature(x = "merMod")}: ... }
    \item{update}{\code{signature(object = "merMod")}: ... }
    \item{VarCorr}{\code{signature(x = "merMod")}: ... }
    \item{vcov}{\code{signature(object = "merMod")}: ... }
  }
}
%%\references{}
%%\author{}

\seealso{\code{\link{lmer}}, \code{\link{glmer}}, \code{\link{nlmer}},
  \code{\linkS4class{reModule}}, \code{\linkS4class{feModule}},
  \code{\linkS4class{respModule}}
}
\examples{
showClass("merMod")
}
\keyword{classes}
