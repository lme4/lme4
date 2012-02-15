\docType{class}
\name{merMod-class}
\alias{anova.merMod}
\alias{coef.merMod}
\alias{deviance.merMod}
\alias{fitted.merMod}
\alias{formula.merMod}
\alias{glmerMod-class}
\alias{lmerMod-class}
\alias{logLik.merMod}
\alias{merMod-class}
\alias{model.frame.merMod}
\alias{model.matrix.merMod}
\alias{nlmerMod-class}
\alias{print.merMod}
\alias{print.summary.mer}
\alias{show.merMod}
\alias{show,merMod-method}
\alias{show.summary.mer}
\alias{summary.merMod}
\alias{summary.summary.mer}
\alias{terms.merMod}
\alias{update.merMod}
\alias{vcov.merMod}
\alias{vcov.summary.mer}
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

