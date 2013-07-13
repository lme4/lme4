\name{residuals.merMod}
\alias{residuals.glmResp}
\alias{residuals.lmResp}
\alias{residuals.merMod}
\title{residuals of merMod objects}
\usage{
  \method{residuals}{merMod} (object,
    type = if (isGLMM(object)) "deviance" else "response",
    scaled = FALSE, ...)

  \method{residuals}{lmResp} (object,
    type = c("working", "response", "deviance", "pearson", "partial"),
    ...)

  \method{residuals}{glmResp} (object,
    type = c("deviance", "pearson", "working", "response", "partial"),
    ...)
}
\arguments{
  \item{object}{a fitted [g]lmer (\code{merMod}) object}

  \item{type}{type of residuals}

  \item{scaled}{scale residuals by residual standard
  deviation (=scale parameter)?}

  \item{\dots}{additional arguments (ignored: for method
  compatibility)}
}
\description{
  residuals of merMod objects
}
\details{
  \itemize{ \item The default residual type varies between
  \code{lmerMod} and \code{glmerMod} objects: they try to
  mimic \code{\link{residuals.lm}} and
  \code{\link{residuals.glm}} respectively. In particular,
  the default \code{type} is \code{"response"}, i.e.
  (observed-fitted) for \code{lmerMod} objects vs.
  \code{"deviance"} for \code{glmerMod} objects.
  \code{type="partial"} is not yet implemented for either
  type. \item Note that the meaning of \code{"pearson"}
  residuals differs between \code{\link{residuals.lm}} and
  \code{\link{residuals.lme}}.  The former returns values
  scaled by the square root of user-specified weights (if
  any), but \emph{not} by the residual standard deviation,
  while the latter returns values scaled by the estimated
  standard deviation (which will include the effects of any
  variance structure specified in the \code{weights}
  argument).  To replicate \code{lme} behaviour, use
  \code{type="pearson"}, \code{scaled=TRUE}. }
}

