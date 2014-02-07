\name{lme4-package}
\alias{lme4}
\alias{lme4-package}
\docType{package}
\title{Linear, generalized linear, and nonlinear mixed models}
\description{
  \code{lme4} provides functions for fitting and analyzing
  mixed models: linear (\code{\link{lmer}}), generalized linear
  (\code{\link{glmer}}) and nonlinear (\code{\link{nlmer}}.)
}
\section{Differences between \pkg{nlme} and \pkg{lme4}}{
  \pkg{lme4} covers approximately the same ground as the earlier
  \pkg{nlme} package.  The most important differences are:
  \itemize{
    \item \pkg{lme4} uses modern, efficient linear algebra methods
    as implemented in the \code{Eigen} package, and uses reference
    classes to avoid undue copying of large objects; it is therefore likely
    to be faster and more memory-efficient than \pkg{nlme}.
    \item \pkg{lme4} includes generalized linear mixed model (GLMM)
    capabilities, via the \code{\link{glmer}} function.
    \item \pkg{lme4} does \emph{not} currently implement \pkg{nlme}'s
    features for modeling heteroscedasticity and
    correlation of residuals.
    \item \pkg{lme4} does not currently offer the same flexibility
    as \pkg{nlme} for composing complex variance-covariance structures,
    but it does implement crossed random effects in a way that
    is both easier for the user and much faster.
    \item \pkg{lme4} offers built-in facilities for likelihood
    profiling and parametric bootstrapping.
    \item \pkg{lme4} is designed to be more modular than \pkg{nlme},
    making it easier for downstream package developers and end-users
    to re-use its components for extensions of the basic mixed model
    framework. It also allows more flexibility for specifying
    different functions for optimizing over the random-effects
    variance-covariance parameters.
    \item \pkg{lme4} is not (yet) as well-documented as \pkg{nlme}.
  }
}
\section{Differences between current (1.0.+) and previous versions of
  \pkg{lme4}}{
  \itemize{
  \item \code{[gn]lmer} now produces objects of class \code{\linkS4class{merMod}}
  rather than class \code{mer} as before
  \item the new version uses a combination of S3 and reference classes
  (see \code{\link{ReferenceClasses}}, \code{\link{merPredD-class}}, and
  \code{\link{lmResp-class}}) as well as S4 classes; partly for this reason
  it is more interoperable with \pkg{nlme}
  \item The internal structure of [gn]lmer is now more modular, allowing
  finer control of the different steps of argument checking; construction
  of design matrices and data structures; parameter estimation; and construction
  of the final \code{merMod} object (see \code{\link{modular}})
  \item profiling and parametric bootstrapping are new in
  the current version
  \item the new version of \pkg{lme4} does \emph{not} provide
  an \code{mcmcsamp} (post-hoc MCMC sampling) method, because this
  was deemed to be unreliable.  Alternatives for computing p-values
  include parametric bootstrapping (\code{\link{bootMer}}) or methods
  implemented in the \pkg{pbkrtest} package and leveraged by the
  \pkg{lmerTest} package and the \code{Anova} function in the \pkg{car} package
  (see \code{\link{pvalues}} for more details).
}
}
\section{Caveats and trouble-shooting}{
\itemize{
    \item Some users who have previously installed versions of the
    RcppEigen and minqa packages may encounter segmentation faults (!!);
    the solution is to make sure to re-install these packages before
    installing \pkg{lme4}.  (Because the problem is not with the
    explicit \emph{version} of the packages, but with running
    packages that were built with different versions of \pkg{Rcpp}
    in conjunction with each other, simply making sure you have
    the latest version, or using \code{update.packages}, will
    not necessarily solve the problem; you must actually re-install
    the packages. The problem is most likely with \pkg{minqa}.)
  }

}
