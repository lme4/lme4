\name{merMod-class}
\title{Class "merMod" of Fitted Mixed-Effect Models}
\Rdversion{1.1}
\docType{class}
\alias{merMod-class}
\alias{lmerMod-class}
\alias{glmerMod-class}
\alias{nlmerMod-class}
%% FIXME --- currently only sigma() and show() ...
\alias{anova,merMod-method}
\alias{coef,merMod-method}
%\alias{devcomp,merMod-method}
\alias{deviance,merMod-method}
\alias{fitted,merMod-method}
%\alias{fixef,merMod-method}
\alias{formula,merMod-method}
\alias{logLik,merMod-method}
\alias{model.frame,merMod-method}
\alias{model.matrix,merMod-method}
\alias{print,merMod-method}
%\alias{profile,merMod-method}
%\alias{ranef,merMod-method}
\alias{show,merMod-method}
%\alias{sigma,merMod-method}-> ./sigma.Rd
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
  A mixed-effects model represented as a \code{\linkS4class{merPredD}} object
  and a response module of a class that inherits from class
  \code{\linkS4class{lmResp}}.  A model with a
  \code{\linkS4class{lmerResp}} response has class \code{lmerMod};
  a \code{\linkS4class{glmResp}} response has class \code{glmerMod};
  and a \code{\linkS4class{nlsResp}} response has class \code{nlmerMod}.
}
\section{Objects from the Class}{
  Objects are created by calls to \code{\link{lmer}},
  \code{\link{glmer}} or \code{\link{nlmer}}.
}
\section{Slots}{
  \describe{
    \item{\code{Gp}:}{The group pointers - an increasing vector of
      non-negative integers that must start with 0.}
    \item{\code{call}:}{The matched call to the function that generated
      the object.}
    \item{\code{frame}:}{The model frame used to evaluate the model matrices.}
    \item{\code{flist}:}{Object of class \code{"list"} ~~ }
    \item{\code{cnms}:}{Object of class \code{"list"} ~~ }
    \item{\code{lower}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{theta}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{beta}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{u}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{devcomp}:}{Object of class \code{"list"} ~~ }
    \item{\code{pp}:}{A prediction module object of class
      \code{"\linkS4class{merPredD}"}.}
    \item{\code{devcomp}:}{A list with two components:
      \dQuote{\code{cmp}}, a numeric vector of components used to
      evaluate the deviance, and \dQuote{\code{dims}}, an integer vector
      of dimensions and flags.}
    \item{\code{re}:}{An \code{\linkS4class{lmerResp}},
      \code{\linkS4class{glmResp}} or \code{\linkS4class{nlsResp}}
      object containing information on the response, the fitted values,
      the residuals and various forms of weights.  According to the
      class of this object the model is a linear mixed model
      (\code{\link{lmer}}), a generalized linear mixed model
      (\code{\link{glmer}}) or a nonlinear mixed model
      (\code{\link{nlmer}}).}
  }
}
\section{Methods}{
  \describe{
%% FIXME: all these are *not* there as S4 methods, apart from  show() and sigma()
    \item{anova}{\code{signature(object = "merMod")}: ... }
    \item{coef}{\code{signature(object = "merMod")}: ... }
    \item{devcomp}{\code{signature(x = "merMod")}: ... }
    \item{deviance}{\code{signature(object = "merMod")}: ... }
    \item{fixef}{\code{signature(object = "merMod")}: ... }
    \item{formula}{\code{signature(x = "merMod")}: ... }
    \item{getME}{\code{signature(x = "merMod")}: ... }
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
  \code{\linkS4class{merPredD}}, \code{\linkS4class{lmerResp}},
  \code{\linkS4class{glmResp}}, \code{\linkS4class{nlsResp}}
}
\examples{
showClass("merMod")
}
\keyword{classes}
