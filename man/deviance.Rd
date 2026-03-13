\name{deviance}
\alias{deviance}
\alias{deviance.merMod}

\title{Deviance of Fitted Mixed-Effect Models}

\description{
  Calculates the deviance of a model of class merMod. Some important details
  regarding the calculation of the deviance is found below.
}

\usage{
\S3method{deviance}{merMod}(object, REML = NULL, \dots)
}

\arguments{
  \item{object}{an \R object of class \code{\linkS4class{merMod}}, i.e.,
    as resulting from \code{\link{lmer}()}, or \code{\link{glmer}()},
    etc.}
   \item{REML}{Logical. If \code{TRUE}, return the restricted log-likelihood
    rather than the log-likelihood.  If \code{NULL} (the default),
    set \code{REML} to \code{isREML(object)} (see \code{\link{isREML}}).}
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
  involves integrating out these coefficients.  The following table
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
  
  Note that \code{summary()} reports the unconditional absolute deviance,
  \code{-2*logLik(object)}.

  For more information about this topic see the \code{misc/logLikGLMM}
  directory in the package source.
}

