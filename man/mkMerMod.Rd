\name{mkMerMod}
\title{Create a 'merMod' Object}
\alias{mkMerMod}
\description{
  Create an object of (a subclass of) class \code{\linkS4class{merMod}}
  from the environment of the objective function and the value returned
  by the optimizer.
}
\usage{
  mkMerMod(rho, opt, reTrms, fr, mc, lme4conv = NULL)
}
\arguments{
  \item{rho}{the environment of the objective function}
  \item{opt}{the optimization result returned by the optimizer
    (a \code{\link{list}}: see \code{\link{lmerControl}} for required elements)}
  \item{reTrms}{random effects structure from the calling function
       (see \code{\link{mkReTrms}} for required elements)}
  \item{fr}{model frame (see \code{\link{model.frame}})}
  \item{mc}{matched call from the calling function}
  \item{lme4conv}{lme4-specific convergence information
    (results of \code{checkConv})}
}
\value{
  an object from a class that inherits from \code{\linkS4class{merMod}}.
}
%% FIXME: TODO
%% \examples{
%% ## 1) An "empty"   merMod object :
%% ## 2) A  "lmer()-like" merMod object, using our "modular" approach instead of lmer()
%% }
