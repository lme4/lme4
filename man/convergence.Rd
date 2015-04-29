\name{convergence}
\alias{convergence}
\title{Assessing convergence for fitted models}
\description{
  The \code{lme4} package uses various implementations
  of general-purpose nonlinear optimizers (e.g. Nelder-Mead or Powell's
  BOBYQA method) to estimate the variance-covariance matrices of
  the random effects.  Assessing reliably whether such
  algorithms have converged is difficult.  For example, evaluating the
  \href{http://en.wikipedia.org/wiki/Karush%E2%80%93Kuhn%E2%80%93Tucker_conditions}{Karush-Kuhn-Tucker conditions} (convergence criteria which in the simplest case of non-constrained optimization reduce to showing that the gradient is zero and the Hessian is positive definite) is challenging because of the difficulty of evaluating the gradient and Hessian.
}
\seealso{\code{\link{lmerControl}}}
