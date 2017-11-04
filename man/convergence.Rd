\name{convergence}
\alias{convergence}
\title{Assessing Convergence for Fitted Models}

\description{
  The \pkg{lme4} package uses general-purpose nonlinear optimizers
  (e.g. Nelder-Mead or Powell's BOBYQA method) to estimate the
  variance-covariance matrices of the random effects.  Assessing
  reliably whether such algorithms have converged is difficult.  For
  example, evaluating the
  \href{http://en.wikipedia.org/wiki/Karush\%E2\%80\%93Kuhn\%E2\%80\%93Tucker_conditions}{%
    Karush-Kuhn-Tucker conditions} (convergence criteria which in the
  simplest case of non-constrained optimization reduce to showing that
  the gradient is zero and the Hessian is positive definite) is
  challenging because of the difficulty of evaluating the gradient and
  Hessian.

  We (the \code{lme4} authors and maintainers) are still in the process
  of finding the best strategies for testing convergence.  Some of the
  relevant issues are
  \itemize{
    \item the gradient and Hessian are the basic ingredients of
    KKT-style testing, but when they have to be estimated by finite
    differences (as in the case of \code{lme4}; direct computation
    of derivatives based on analytic expressions may eventually be
    available for some special classes, but we have not yet
    implemented them) they may not be sufficiently accurate for
    reliable convergence testing.

    \item The Hessian computation in particular represents
    a difficult tradeoff between computational expense and
    accuracy.  At present the Hessian computations used
    for convergence checking (and for estimating standard errors
    of fixed-effect parameters for GLMMs) follow the \pkg{ordinal} package
    in using a naive but computationally cheap centered finite difference
    computation (with a fixed step size of \eqn{10^{-4}}{1e-4}).  A more
    reliable but more expensive approach is to use
    \href{http://en.wikipedia.org/wiki/Richardson_extrapolation}{Richardson extrapolation},
    as implemented in the \pkg{numDeriv} package.

    \item it is important to scale the estimated gradient at
    the estimate appropriately; two reasonable approaches are
    \enumerate{
      \item don't scale random-effects (Cholesky) gradients, since these
      are essentially already unitless (for LMMs they are scaled
      relative to the residual variance; for GLMMs they are scaled
      relative to the sampling variance of the conditional distribution);
      for GLMMs, scale fixed-effect gradients by the standard deviations
      of the corresponding input variable, or
      \item scale gradients by the inverse Cholesky factor of the
      Hessian, equivalent to scaling by the estimated Wald standard error
      of the estimated parameters.  The latter approach is used in the
      current version of \pkg{lme4}; it has the disadvantage that it
      requires us to estimate the Hessian (although the Hessian is
      required \href{https://github.com/lme4/lme4/issues/47}{for
	reliable estimation of the fixed-effect standard errors for GLMMs}
      in any case).
    }
    \item Exploratory analyses suggest that (1) the naive estimation
    of the Hessian may fail for large data sets (number of observations
    greater than approximately
    \eqn{10^{5}}{1e5}); (2) the magnitude of the scaled
    gradient increases with sample size, so that warnings will occur
    even for apparently well-behaved fits with large data sets.
  }
  If you do see convergence warnings, and want to trouble-shoot/double-check
  the results, the following steps are recommended (examples are given below):
  \itemize{
    \item double-check the model specification and the data
    for mistakes
    \item center and scale continuous predictor variables (e.g. with
    \code{\link{scale}})
    \item check for singularity: if any of the diagonal elements of the
    Cholesky factor are zero or very small, the convergence testing methods may be
    inappropriate (see examples)
    \item double-check the Hessian calculation with the more expensive
    Richardson extrapolation method (see examples)
    \item restart the fit from the apparent optimum, or from a point
    perturbed slightly away from the optimum
    \item try all available optimizers (e.g. several different implementations
    of BOBYQA and Nelder-Mead, L-BFGS-B from \code{optim}, \code{nlminb},
    \dots) via the \code{allFit()} function, see \sQuote{5.} in the
    examples.  While this will of course be slow for large fits, we consider
    it the gold standard; if all optimizers converge to values that
    are practically equivalent, then we would consider the convergence
    warnings to be false positives.
  }
  To quote Douglas Adams,
  \href{http://en.wikipedia.org/wiki/So_Long,_and_Thanks_for_All_the_Fish}{we
    apologize for the inconvenience}.
}% end{description}

\examples{
fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

## 1. center and scale predictors:
ss.CS <- transform(sleepstudy, Days=scale(Days))
fm1.CS <- update(fm1, data=ss.CS)

## 2. check singularity
diag.vals <- getME(fm1,"theta")[getME(fm1,"lower") == 0]
any(diag.vals < 1e-6) # FALSE

## 3. recompute gradient and Hessian with Richardson extrapolation
devfun <- update(fm1, devFunOnly=TRUE)
if (isLMM(fm1)) {
    pars <- getME(fm1,"theta")
} else {
    ## GLMM: requires both random and fixed parameters
    pars <- getME(fm1, c("theta","fixef"))
}
if (require("numDeriv")) {
    cat("hess:\n"); print(hess <- hessian(devfun, unlist(pars)))
    cat("grad:\n"); print(grad <- grad(devfun, unlist(pars)))
    cat("scaled gradient:\n")
    print(scgrad <- solve(chol(hess), grad))
}
## compare with internal calculations:
fm1@optinfo$derivs

## 4. restart the fit from the original value (or
## a slightly perturbed value):
fm1.restart <- update(fm1, start=pars)

## 5. try all available optimizers

  source(system.file("utils", "allFit.R", package="lme4"))
  fm1.all <- allFit(fm1)
  ss <- summary(fm1.all)
  ss$ fixef               ## extract fixed effects
  ss$ llik                ## log-likelihoods
  ss$ sdcor               ## SDs and correlations
  ss$ theta               ## Cholesky factors
  ss$ which.OK            ## which fits worked

}

\seealso{\code{\link{lmerControl}}}
