\name{glmerLaplaceHandle}
\alias{glmerLaplaceHandle}
\title{Handle for \code{glmerLaplace}}
\usage{
glmerLaplaceHandle(pp, resp, nAGQ, tol, maxit, verbose,
                   dispProfile = TRUE, maxPhiIter = 100L, qEff = NA_real_)
}
\arguments{
  \item{pp}{\code{\link{merPredD}} object}
  \item{resp}{\code{\link{lmResp}} object}
  \item{nAGQ}{see \code{\link{glmer}}}
  \item{tol}{tolerance}
  \item{maxit}{maximum number of pwrss iterations}
  \item{verbose}{display optimizer progress}
  \item{dispProfile}{for families with a free dispersion parameter
    (e.g. Gamma, gaussian, inverse.gaussian), profile the dispersion
    parameter via a nested fixed-point loop within PIRLS (see
    \code{\link{glmerControl}}'s \code{disp_method})}
  \item{maxPhiIter}{maximum number of iterations of the nested
    dispersion-profiling loop (only used when \code{dispProfile=TRUE})}
  \item{qEff}{degrees-of-freedom correction to the moment estimator's
    denominator (see \code{\link{glmerControl}}'s
    \code{disp_dof_correction}), or \code{NA_real_} for none}
}
\description{
  Handle for calling the \code{glmerLaplace} \code{C++} function.  Not
  intended for routine use.
}
\value{
  Value of the objective function
}
\keyword{utilities}
