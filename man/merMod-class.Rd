\docType{class}
\name{merMod-class}
\alias{anova,merMod-method}
\alias{coef,merMod-method}
\alias{deviance,merMod-method}
\alias{fitted,merMod-method}
\alias{formula,merMod-method}
\alias{glmerMod-class}
\alias{lmerMod-class}
\alias{logLik,merMod-method}
\alias{merMod-class}
\alias{model.frame,merMod-method}
\alias{model.matrix,merMod-method}
\alias{nlmerMod-class}
\alias{print,merMod-method}
\alias{print,summary.mer-method}
\alias{show,merMod-method}
\alias{show,summary.mer-method}
\alias{simulate,merMod-method}
\alias{summary,merMod-method}
\alias{summary,summary.mer-method}
\alias{terms,merMod-method}
\alias{update,merMod-method}
\alias{VarCorr,merMod-method}
\alias{vcov,merMod-method}
\alias{vcov,summary.mer-method}
\title{Class "merMod" of Fitted Mixed-Effect Models}
\description{
  A mixed-effects model represented as a
  \code{\linkS4class{merPredD}} object and a response
  module of a class that inherits from class
  \code{\linkS4class{lmResp}}.  A model with a
  \code{\linkS4class{lmerResp}} response has class
  \code{lmerMod}; a \code{\linkS4class{glmResp}} response
  has class \code{glmerMod}; and a
  \code{\linkS4class{nlsResp}} response has class
  \code{nlmerMod}.
}
\section{Objects from the Class}{
  Objects are created by calls to \code{\link{lmer}},
  \code{\link{glmer}} or \code{\link{nlmer}}.
}
\examples{
showClass("merMod")
}
\seealso{
  \code{\link{lmer}}, \code{\link{glmer}},
  \code{\link{nlmer}}, \code{\linkS4class{merPredD}},
  \code{\linkS4class{lmerResp}},
  \code{\linkS4class{glmResp}},
  \code{\linkS4class{nlsResp}}
}
\keyword{classes}

